# Snakefile

configfile: "inputs/config.yaml"

EXPERIMENTS = list(config["experiments"].keys())
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
        expand("results/{exp}/osdmap.bin", exp=EXPERIMENTS),
        expand("results/{exp}/osdmap.txt", exp=EXPERIMENTS),
        expand("results/{exp}/crushmap.bin", exp=EXPERIMENTS),
        expand("results/{exp}/crushmap.txt", exp=EXPERIMENTS),
        expand("results/{exp}/crushmap.json", exp=EXPERIMENTS),


rule build_maps:
    input:
        unpack(get_inputs)
    output:
        osdmap_bin = "results/{exp}/osdmap.bin",
        osdmap_txt = "results/{exp}/osdmap.txt",
        crushmap_bin = "results/{exp}/crushmap.bin",
        crushmap_txt = "results/{exp}/crushmap.txt",
        crushmap_json = "results/{exp}/crushmap.json",
        pgdump = "results/{exp}/pgdump.txt",
    params:
        container_runtime = config["container_runtime"],
        mon_bootstrap_script = config["mon_bootstrap_script"],
        ceph_image = config["ceph_image"],
        out_dir = "results/{exp}",
    script:
        "tools/gen_osdmap.py"
