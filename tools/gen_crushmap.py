from jinja2 import Environment, FileSystemLoader
import os

template_dir = os.path.dirname(snakemake.input.template)
template_file = os.path.basename(snakemake.input.template)
environment = Environment(loader=FileSystemLoader(template_dir))
template = environment.get_template(template_file)

class Node:
    """A node class for a tree structure"""

    def __init__(self, prefix: str, weight: float, children: list['Node'] | None = None):
        self.weight_ = weight
        self.children = children or []
        if weight > 0:
            self.id = get_osd_id()
        else:
            self.id = get_bucket_id()
            self.hdd_id = get_bucket_id()
        self.name = f"{prefix}{self.id}"
    
    def _calculate_weight(self):
        total = self.weight_
        for child in self.children:
            total += child._calculate_weight()
        return total

    def weight(self):
        """Calculates total weight of self + children"""
        raw = self._calculate_weight()
        return f"{raw:.5f}"
    
    def add_child(self, child):
        self.children.append(child)
    
    def get_osds(self):
        """resursive get all childrens children"""
        all_osds = []
        for child in self.children:
            if child.weight_ > 0:
                all_osds.append(child)
            else:
                all_osds.extend(child.get_osds())
        return all_osds


# first available bucket id
bucket_id = -1
def get_bucket_id():
    global bucket_id
    previous_id = bucket_id
    bucket_id -= 1
    return previous_id

# first available osd id
osd_id = 0
def get_osd_id():
    global osd_id
    previous_id = osd_id
    osd_id += 1
    return previous_id


osds_per_host_num = snakemake.params.osds_per_host_num
hosts_per_rack_table = snakemake.params.num_hosts_per_rack_table
racks_per_dc_table = snakemake.params.num_racks_per_dc_table
osd_weight_dc_table = snakemake.params.osd_weight_dc_table


def gen_crushmap():
    root = Node("root", 0)
    for dc_i in range(len(racks_per_dc_table)):
        dc_str = f"dc{dc_i}"
        dc = Node("dc", 0)
        osd_weight = osd_weight_dc_table[dc_i]

        for rack_i in range(racks_per_dc_table[dc_i]):
            rack_str = f"rack{rack_i}"
            rack = Node("rack", 0)

            for host_i in range(hosts_per_rack_table[rack_i]):
                host_str = f"host{host_i}"
                host = Node("host", 0)

                for osd_i in range(osds_per_host_num):
                    osd_str = f"osd.{osd_i}"
                    osd = Node("osd.", osd_weight)
                    host.add_child(osd)

                rack.add_child(host)

            dc.add_child(rack)

        root.add_child(dc)


    ctx = {"roots": [root]}

    content = template.render(ctx)

    with open(snakemake.output[0], "w") as f:
        f.write(content)

    print("hello from gen_crushmap.py")

gen_crushmap()
