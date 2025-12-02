# Snakefile

configfile: "inputs/config.yaml"

EXPERIMENTS = list(config["experiments"].keys())
CRUSH_SPEC_PARAMS_IDX = [0,1]
# RESULT_ARTIFACTS = ["osdmap.bin", "osdmap.txt", "crushmap.bin", "crushmap.txt", "crushmap.json"]


# Helpers
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
        osds_per_host_num = lambda w: config["experiments"][w.exp]["crush_spec_params"][int(w.i)]["osds_per_host_num"],
        num_hosts_per_rack_table = lambda w: config["experiments"][w.exp]["crush_spec_params"][int(w.i)]["hosts_per_rack_table"],
        num_racks_per_dc_table = lambda w: config["experiments"][w.exp]["crush_spec_params"][int(w.i)]["racks_per_dc_table"],
        osd_weight_dc_table = lambda w: config["experiments"][w.exp]["crush_spec_params"][int(w.i)]["osd_weight_per_dc_table"],
    script:
        "tools/gen_crushmap.py"
