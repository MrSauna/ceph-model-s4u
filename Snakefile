# Snakefile
import glob

configfile: "inputs/config.yaml"

if "active_experiments" in config and config["active_experiments"]:
    EXPERIMENTS = [e for e in config["active_experiments"] if e in config["experiments"]]
else:
    EXPERIMENTS = list(config["experiments"].keys())
CRUSH_SPEC_PARAMS_IDX = [0,1]
# RESULT_ARTIFACTS = ["osdmap.bin", "osdmap.txt", "crushmap.bin", "crushmap.txt", "crushmap.json"]


# Helpers
def ensure_list(obj):
    if isinstance(obj, list):
        return obj
    return [obj]

def get_recipe_for_exp(wildcards):
    return config["experiments"][wildcards.exp]["recipe"]

def get_inputs(wildcards):
    exp_data = config["experiments"][wildcards.exp]
    inputs = {
        "recipe": exp_data["recipe"]
    }
    
    # Only add the dependency if 'base_map' is defined in config
    if "base_map" in exp_data:
        inputs["base_map"] = exp_data["base_map"]
        
    return inputs


# Rules
rule all:
    input:
        expand("results/{exp}/{i}/osdmap.bin", exp=EXPERIMENTS, i=CRUSH_SPEC_PARAMS_IDX),
        expand("results/{exp}/{i}/osdmap.txt", exp=EXPERIMENTS, i=CRUSH_SPEC_PARAMS_IDX),
        expand("results/{exp}/{i}/crushmap.bin", exp=EXPERIMENTS, i=CRUSH_SPEC_PARAMS_IDX),
        expand("results/{exp}/{i}/crushmap.txt", exp=EXPERIMENTS, i=CRUSH_SPEC_PARAMS_IDX),
        expand("results/{exp}/{i}/crushmap.json", exp=EXPERIMENTS, i=CRUSH_SPEC_PARAMS_IDX),

        expand("results/{exp}/net_metrics.csv", exp=EXPERIMENTS),
        expand("results/{exp}/mon_metrics.csv", exp=EXPERIMENTS),
        expand("results/{exp}/client_metrics.csv", exp=EXPERIMENTS),
        expand("results/{exp}/topology.json", exp=EXPERIMENTS),

        expand("results/{exp}/figures/pgmap_plot.svg", exp=EXPERIMENTS),
        expand("results/{exp}/figures/client_throughput.svg", exp=EXPERIMENTS),
        expand("results/{exp}/figures/network_topology.svg", exp=EXPERIMENTS),


rule viz_pgmap:
    input:
        net_metrics = "results/{exp}/net_metrics.csv",
        mon_metrics = "results/{exp}/mon_metrics.csv",
        client_metrics = "results/{exp}/client_metrics.csv",
    output:
        client_throughput = "results/{exp}/figures/client_throughput.svg",
        pgmap_plot = "results/{exp}/figures/pgmap_plot.svg",
        network_topology = "results/{exp}/figures/network_topology.svg",
    params:
        output_dir = "results/{exp}/figures",
    script: 
        "tools/visualize.py"

rule run_sim:
    input:
        binary = "build_release/ceph-sim",
        pgdump1 = "results/{exp}/0/pgdump.txt",
        pgdump2 = "results/{exp}/1/pgdump.txt",
    output:
        net_metrics = "results/{exp}/net_metrics.csv",
        mon_metrics = "results/{exp}/mon_metrics.csv",
        client_metrics = "results/{exp}/client_metrics.csv",
        topology = "results/{exp}/topology.json",
        stdout = "results/{exp}/stdout.txt",
        stderr = "results/{exp}/stderr.txt",
    params:
        clients_num = lambda w: config["experiments"][w.exp]["clients_num"],
        profile = lambda w: config["experiments"][w.exp].get("profile", "balanced"),
        client_read_queue_depth = lambda w: config["experiments"][w.exp].get("client_read_queue_depth", 1),
        client_write_queue_depth = lambda w: config["experiments"][w.exp].get("client_write_queue_depth", 1),
        dc_shape_csv = lambda w: ",".join(
            ensure_list(config["experiments"][w.exp]["crush_spec_params"][-1].get("dc_shape", ["2:2:2"]))
        ),
        dc_speed_csv = lambda w: ",".join(
            ensure_list(config["experiments"][w.exp]["crush_spec_params"][-1].get("dc_speed", ["10:10:10"]))
        ),
    shell:
        """
        {input.binary} \
            "--dc-shape={params.dc_shape_csv}" \
            "--dc-speed={params.dc_speed_csv}" \
            "--pgdump={input.pgdump1}" \
            "--pgdump={input.pgdump2}" \
            "--profile={params.profile}" \
            "--pg-objects=1000" \
            "--object-size=4194304" \
            "--pool-id=1" \
            "--disk-read-bandwidth=125829120" \
            "--disk-write-bandwidth=83886080" \
            "--dc-clients={params.clients_num}" \
            "--client-read-queue-depth={params.client_read_queue_depth}" \
            "--client-write-queue-depth={params.client_write_queue_depth}" \
            "--output-dir=results/{wildcards.exp}" \
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
        crushmap_txt = "results/{exp}/{i}/crushmap.txt",
        mon_bootstrap_script = "inputs/mon_bootstrap.sh",
    output:
        osdmap_bin = "results/{exp}/{i}/osdmap.bin",
        osdmap_txt = "results/{exp}/{i}/osdmap.txt",
        crushmap_bin = "results/{exp}/{i}/crushmap.bin",
        crushmap_json = "results/{exp}/{i}/crushmap.json",
        pgdump = "results/{exp}/{i}/pgdump.txt",
    params:
        container_runtime = config["container_runtime"],
        mon_bootstrap_script = config["mon_bootstrap_script"],
        ceph_image = config["ceph_image"],
        out_dir = "results/{exp}/{i}",
        osd_num = lambda w: config["experiments"][w.exp]["osd_num"],
        pg_num = lambda w: config["experiments"][w.exp]["pg_num"],
    script:
        "tools/gen_osdmap.py"


rule gen_crushmap:
    input:
        template = lambda w: config["experiments"][w.exp]["template"],
        osds_template = "inputs/templates/osds.j2",
        buckets_template = "inputs/templates/buckets.j2",
    output:
        crushmap_txt = "results/{exp}/{i}/crushmap.txt",
    params:
        dc_shape = lambda w: config["experiments"][w.exp]["crush_spec_params"][int(w.i)]["dc_shape"],
        dc_weight = lambda w: config["experiments"][w.exp]["crush_spec_params"][int(w.i)]["dc_weight"],
        osd_num = lambda w: config["experiments"][w.exp]["osd_num"],
    script:
        "tools/gen_crushmap.py"
