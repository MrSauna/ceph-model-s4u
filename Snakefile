# Snakefile
import glob
import subprocess


configfile: "inputs/config.yaml"

if "active_experiments" in config and config["active_experiments"]:
    CASES = [e for e in config["active_experiments"] if e in config["cases"]]
else:
    CASES = list(config["cases"].keys())
CRUSH_SPEC_PARAMS_IDX = [0,1]

SCENARIOS = list(config["scenarios"].keys())

# get git hash
try:
    GIT_HASH = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], 
        cwd=workflow.basedir  # Ensures you run git in the workflow directory
    ).decode("utf-8").strip()
    
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=workflow.basedir):
        GIT_HASH += "-dirty"
except subprocess.CalledProcessError:
    GIT_HASH = "unknown"


# Helpers
def ensure_list(obj):
    if isinstance(obj, list):
        return obj
    return [obj]


def get_recipe_for_casee(wildcards):
    return config["cases"][wildcards.case]["recipe"]


def get_case_inputs(wildcards):
    exp_data = config["cases"][wildcards.case]
    inputs = {
        "recipe": exp_data["recipe"]
    }
    
    # Only add the dependency if 'base_map' is defined in config
    if "base_map" in exp_data:
        inputs["base_map"] = exp_data["base_map"]
        
    return inputs


def get_scenario_inputs(wildcards):
    cases = config["scenarios"][wildcards.scenario]

    return expand("results/{case}/done", case=cases)


wildcard_constraints:
    case = "|".join(CASES),
    scenario = "|".join(SCENARIOS)

# Rules
rule all:
    input:
        expand("results/{case}/done", case=CASES),
        expand("results/{scenario}/done", scenario=SCENARIOS),

rule scenario:
    output:
        done = "results/{scenario}/done",
    input:
        comparison = "results/{scenario}/figures/client_throughput.svg",
    shell:
        """
        touch {output.done}
        """

rule scenario_viz:
    input:
        client_metrics = lambda w: expand("results/{case}/client_metrics.csv", case=config["scenarios"][w.scenario]),
        net_metrics = lambda w: expand("results/{case}/net_metrics.csv", case=config["scenarios"][w.scenario]),
        mon_metrics = lambda w: expand("results/{case}/mon_metrics.csv", case=config["scenarios"][w.scenario]),
        script = "tools/compare_results.py",
    output:
        comparison = "results/{scenario}/figures/client_throughput.svg",
    params:
        git_hash = GIT_HASH,
        output_dir = "results/{scenario}/figures",
    script:
        "tools/compare_results.py"


rule case:
    output:
        done = "results/{case}/done",
    input:
        client_throughput = "results/{case}/figures/client_throughput.svg",
        pgmap_plot = "results/{case}/figures/pgmap_plot.svg",
        network_topology = "results/{case}/figures/network_topology.svg",
        metadata = "results/{case}/metadata.json",
    shell:
        """
        touch {output.done}
        """


rule metadata:
    output:
        metadata = "results/{case}/metadata.json",
    params:
        case = lambda w: w.case,
        git_hash = GIT_HASH,
    shell:
        """
        echo '{{"case": "{params.case}", "git_hash": "{params.git_hash}"}}' > {output.metadata}
        """


rule viz_pgmap:
    input:
        net_metrics = "results/{case}/net_metrics.csv",
        mon_metrics = "results/{case}/mon_metrics.csv",
        client_metrics = "results/{case}/client_metrics.csv",
    output:
        client_throughput = "results/{case}/figures/client_throughput.svg",
        pgmap_plot = "results/{case}/figures/pgmap_plot.svg",
        network_topology = "results/{case}/figures/network_topology.svg",

    params:
        output_dir = "results/{case}/figures",
        git_hash = GIT_HASH,
    script: 
        "tools/visualize.py"


rule run_sim:
    input:
        binary = "build_release/ceph-sim",
        pgdump1 = "results/{case}/0/pgdump.txt",
        pgdump2 = "results/{case}/1/pgdump.txt",
    output:
        net_metrics = "results/{case}/net_metrics.csv",
        mon_metrics = "results/{case}/mon_metrics.csv",
        client_metrics = "results/{case}/client_metrics.csv",
        topology = "results/{case}/topology.json",
        stdout = "results/{case}/stdout.txt",
        stderr = "results/{case}/stderr.txt",
    params:
        clients_num = lambda w: config["cases"][w.case]["clients_num"],
        start_up_delay = lambda w: config["cases"][w.case].get("start_up_delay", 0),
        shut_down_delay = lambda w: config["cases"][w.case].get("shut_down_delay", 0),
        profile = lambda w: config["cases"][w.case].get("profile", "balanced"),
        client_read_queue_depth = lambda w: config["cases"][w.case].get("client_read_queue_depth", 1),
        client_write_queue_depth = lambda w: config["cases"][w.case].get("client_write_queue_depth", 1),
        dc_shape_csv = lambda w: ",".join(
            ensure_list(config["cases"][w.case]["crush_spec_params"][-1].get("dc_shape", ["2:2:2"]))
        ),
        dc_speed_csv = lambda w: ",".join(
            ensure_list(config["cases"][w.case]["crush_spec_params"][-1].get("dc_speed", ["10:10:10"]))
        ),
        pg_objects = lambda w: config["cases"][w.case].get("pg_objects", 1000),
        object_size = lambda w: config["cases"][w.case].get("object_size", 4194304),
    shell:
        """
        {input.binary} \
            "--dc-shape={params.dc_shape_csv}" \
            "--dc-speed={params.dc_speed_csv}" \
            "--pgdump={input.pgdump1}" \
            "--pgdump={input.pgdump2}" \
            "--start-up-delay={params.start_up_delay}" \
            "--shut-down-delay={params.shut_down_delay}" \
            "--profile={params.profile}" \
            "--pg-objects={params.pg_objects}" \
            "--object-size={params.object_size}" \
            "--pool-id=1" \
            "--disk-write-bandwidth=209715200" \
            "--dc-clients={params.clients_num}" \
            "--client-read-queue-depth={params.client_read_queue_depth}" \
            "--client-write-queue-depth={params.client_write_queue_depth}" \
            "--output-dir=results/{wildcards.case}" \
            > {output.stdout} 2> {output.stderr}
        """


SOURCES = glob.glob("src/*.cpp") + glob.glob("src/*.hpp")

rule compile_sim:
    input:
        src = SOURCES,
        cmakelists = "CMakeLists.txt",
    output:
        binary = "build_release/ceph-sim",
    params:
        build_dir = "build_release",
    shell:
        """
        mkdir -p {params.build_dir}
        cmake -S . -B {params.build_dir} -DCMAKE_BUILD_TYPE=Release
        cmake --build {params.build_dir} --target ceph-sim -- -j2
        """


rule build_maps:
    input:
        crushmap_txt = "results/{case}/{i}/crushmap.txt",
        mon_bootstrap_script = "inputs/mon_bootstrap.sh",
    output:
        osdmap_bin = "results/{case}/{i}/osdmap.bin",
        osdmap_txt = "results/{case}/{i}/osdmap.txt",
        crushmap_bin = "results/{case}/{i}/crushmap.bin",
        crushmap_json = "results/{case}/{i}/crushmap.json",
        pgdump = "results/{case}/{i}/pgdump.txt",
    params:
        container_runtime = config["container_runtime"],
        mon_bootstrap_script = config["mon_bootstrap_script"],
        ceph_image = config["ceph_image"],
        out_dir = "results/{case}/{i}",
        osd_num = lambda w: config["cases"][w.case]["osd_num"],
        pg_num = lambda w: config["cases"][w.case]["pg_num"],
    script:
        "tools/gen_osdmap.py"


rule gen_crushmap:
    input:
        template = lambda w: config["cases"][w.case]["template"],
        osds_template = "inputs/templates/osds.j2",
        buckets_template = "inputs/templates/buckets.j2",
    output:
        crushmap_txt = "results/{case}/{i}/crushmap.txt",
    params:
        dc_shape = lambda w: config["cases"][w.case]["crush_spec_params"][int(w.i)]["dc_shape"],
        dc_weight = lambda w: config["cases"][w.case]["crush_spec_params"][int(w.i)]["dc_weight"],
        osd_num = lambda w: config["cases"][w.case]["osd_num"],
    script:
        "tools/gen_crushmap.py"
