import subprocess
import os
import uuid
import time
import sys


def setup_monitor(*, container_runtime, ceph_image, c_name, volumes, mon_bootstrap_script):

    abs_bootstrap = os.path.abspath(mon_bootstrap_script)
    volumes.extend([
        "-v", f"{abs_bootstrap}:/mon_bootstrap.sh:z"
    ])

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


    return c_name


def gen_osdmap():
    # Access Parameters from Snakemake
    output_dir = snakemake.params.out_dir
    container_runtime = snakemake.params.container_runtime
    ceph_image = snakemake.params.ceph_image
    mon_bootstrap_script = snakemake.input.mon_bootstrap_script
    crushmap_txt = snakemake.input.crushmap_txt
    osd_num = snakemake.params.osd_num
    pg_num = snakemake.params.pg_num
    
    # Prepare Environment
    os.makedirs(output_dir, exist_ok=True)
    
    abs_crushmap_txt = os.path.abspath(crushmap_txt)
    abs_out = os.path.abspath(output_dir)
    
    # Build volume bindings
    volumes = [
        "-v", f"{abs_out}:/out:z",
        "-v", f"{abs_crushmap_txt}:/cm.txt:z",
    ]
    
    # If we have a base map, mount it to a known location


    c_name = f"ceph_gen_{uuid.uuid4().hex[:8]}"
    try:
        # Start Container
        setup_monitor(container_runtime=container_runtime,
                      ceph_image=ceph_image,
                      c_name=c_name,
                      volumes=volumes,
                      mon_bootstrap_script=mon_bootstrap_script)



        # print("--- Generate osdmap ---")
        # subprocess.run([container_runtime, "exec", c_name, "osdmaptool", "/om", "--createsimple", str(osd_num), "--with-default-pool", "--pg-bits", str(pg_bits), "--clobber"], check=True)
        # subprocess.run([container_runtime, "exec", c_name, "osdmaptool", "/om", "--tree"], check=True)


        print("--- Import crushmap ---")
        subprocess.run([container_runtime, "exec", c_name, "crushtool", "-c", "/cm.txt", "-o", "/cm"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "setcrushmap", "-i", "/cm"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "setmaxosd", str(osd_num)], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "tree"], check=True)

        print("--- Create pool ---")
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "pool", "create", "mypool", str(pg_num)], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "pool", "set", "mypool", "size", "3"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "pool", "set", "mypool", "min_size", "2"], check=True)

        print("--- Exporting Artifacts ---")
        # osdmap
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "getmap", "-o", "/out/osdmap.bin"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "dump", "-o", "/out/osdmap.txt"], check=True)

        # crushmap
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "getcrushmap", "-o", "/out/crushmap.bin"], check=True)
        subprocess.run([container_runtime, "exec", c_name, "ceph", "osd", "crush", "dump", "-o", "/out/crushmap.json"], check=True)

        # pgdump from osdmaptool
        with open(abs_out + "/pgdump.txt", "w") as f:
            subprocess.run([container_runtime, "exec", c_name, "osdmaptool", "/out/osdmap.bin", "--test-map-pgs-dump-all", "--mark-up-in"], stdout=f, check=True)

    finally:
        subprocess.run([container_runtime, "stop", "--time", "1", c_name], stderr=subprocess.DEVNULL)


gen_osdmap()
