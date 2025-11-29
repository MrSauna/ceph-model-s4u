import subprocess
import os
import uuid
import time
import sys


def setup_monitor(container_runtime, ceph_image, c_name, volumes, bootstrap_script, base_osdmap_path=None):

    abs_bootstrap = os.path.abspath(bootstrap_script)
    volumes.extend([
        "-v", f"{abs_bootstrap}:/mon_bootstrap.sh"
    ])

    if base_osdmap_path:
        abs_base = os.path.abspath(base_osdmap_path)
        volumes.extend(["-v", f"{abs_base}:/base_osdmap.bin"])

    print(f"--- Starting Container {c_name} ---")
    subprocess.run([
        container_runtime, "run", "-d", "--rm", "--name", c_name,
        "-e", "MON_IP=127.0.0.1",
        "-e", "CEPH_DAEMON=MON",
        *volumes,
        ceph_image,
        "sleep", "infinity"
    ], check=True)

    # run bootstrap script
    print("--- Bootstrapping Monitor ---")
    subprocess.run([
        container_runtime, "exec", c_name, "bash", "/mon_bootstrap.sh"
    ], check=True)

    # Wait for Health
    print("--- Waiting for Monitor ---")
    for _ in range(30):
        if subprocess.run(
            [container_runtime, "exec", c_name, "ceph", "-s"], 
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        ).returncode == 0:
            break
        time.sleep(1)
    else:
        raise RuntimeError("Ceph Monitor failed to start")

    
    # set base osdmap if provided
    if base_osdmap_path:
        subprocess.run([container_runtime, "exec", c_name, "ceph", "setmap" "-i", "/base_osdmap.bin"], check=True)

    return c_name


def gen_osdmap():
    # Access Parameters from Snakemake
    recipe_path = snakemake.input.recipe
    output_dir = snakemake.params.out_dir
    container_runtime = snakemake.params.container_runtime
    ceph_image = snakemake.params.ceph_image
    mon_bootstrap_script = snakemake.params.mon_bootstrap_script
    
    # Check if optional input exists
    base_osdmap_path = getattr(snakemake.input, "base_map", None)

    # Prepare Environment
    os.makedirs(output_dir, exist_ok=True)
    
    abs_out = os.path.abspath(output_dir)
    abs_recipe = os.path.abspath(recipe_path)
    
    # Build volume bindings
    volumes = [
        "-v", f"{abs_out}:/out",
        "-v", f"{abs_recipe}:/recipe.sh",
    ]
    
    # If we have a base map, mount it to a known location


    c_name = f"ceph_gen_{uuid.uuid4().hex[:8]}"
    try:
        # Start Container
        setup_monitor(container_runtime, ceph_image, c_name, volumes, mon_bootstrap_script)

        # Logic branch: load base map OR start fresh
        if base_osdmap_path:
            print("--- Importing Base Map ---")
            subprocess.run(
                [container_runtime, "exec", c_name, "ceph", "osd", "setmap", "-i", "/base_map.bin"],
                check=True
            )
        
        print(f"--- Running Recipe {recipe_path} ---")
        subprocess.run([container_runtime, "exec", c_name, "bash", "/recipe.sh"], check=True)

        print("--- Exporting Artifacts ---")
        
        # osdmap
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "getmap", "-o", "/out/osdmap.bin"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "dump", "-o", "/out/osdmap.txt"], check=True)

        # crushmap
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "getcrushmap", "-o", "/out/crushmap.bin"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "crushtool", "-d", "/out/crushmap.bin", "-o", "/out/crushmap.txt"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "crush", "dump", "-o", "/out/crushmap.json"], check=True)

        # pgdump from osdmaptool
        with open(abs_out + "/pgdump.txt", "w") as f:
            subprocess.run([container_runtime, "exec", c_name, "osdmaptool", "/out/osdmap.bin", "--test-map-pgs-dump-all"], stdout=f, check=True)

    finally:
        pass
        # subprocess.run([container_runtime, "stop", c_name], stderr=subprocess.DEVNULL)


gen_osdmap()
