from jinja2 import Environment, FileSystemLoader

environment = Environment(loader=FileSystemLoader("templates/"))
template = environment.get_template("pohja.j2")

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
bucket_id = -131
def get_bucket_id():
    global bucket_id
    previous_id = bucket_id
    bucket_id -= 1
    return previous_id

# first available osd id
osd_id = 1227
def get_osd_id():
    global osd_id
    previous_id = osd_id
    osd_id += 1
    return previous_id

region = Node("region", 0)
osd_weights = [22.79800, 20] # 24TB + (4*6.4TB/24) => TiB
for dc_i in range(1):
    dc_str = f"dc{dc_i}"
    dc = Node("dc", 0)
    osd_weight = osd_weights[dc_i]

    for rack_i in range(3):
        rack_str = f"rack{rack_i}"
        rack = Node("rack", 0)

        for host_i in range(2):
            host_str = f"host{host_i}"
            host = Node("host", 0)

            for osd_i in range(24):
                osd_str = f"osd.{osd_i}"
                osd = Node("osd.", osd_weight)
                host.add_child(osd)

            rack.add_child(host)

        dc.add_child(rack)

    region.add_child(dc)


ctx = {"region": region}

content = template.render(ctx)

print(content)
