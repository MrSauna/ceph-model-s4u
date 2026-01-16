from jinja2 import Environment, FileSystemLoader
import os

# --- Snakemake Mocking ---
if "snakemake" not in globals():
    class MockSnakemake:
        def __init__(self):
            self.input = type('Input', (), {})()
            self.input.template = "inputs/templates/cm_simple.j2"
            self.output = ["results/vscode/crushmap.txt"]
            self.params = type('Params', (), {})()
            self.params.dc_shape = ["2:2:2"]
            self.params.dc_weight = ["14.0::"]
            self.params.osd_num = 8
    snakemake = MockSnakemake()

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

# --- Helper Functions (Ported from SimContext.cpp) ---

def split_str(s: str, delim: str) -> list[str]:
    if not s:
        return [""]
    return s.split(delim)

def resolve_val(spec: str, self_id: int, parent_id: int, fallback: float) -> float:
    if not spec:
        return fallback
    
    spec_val = spec
    use_parent = False
    if spec_val.startswith('@'):
        use_parent = True
        spec_val = spec_val[1:]
    
    values = []
    try:
        parts = spec_val.split(',')
        for part in parts:
            if part:
                values.append(float(part))
    except ValueError:
        return fallback

    if not values:
        return fallback
    
    idx = parent_id if use_parent else self_id
    return values[idx % len(values)]

def get_hierarchy_spec(vec: list[str], dc_idx: int, level: int) -> str:
    if dc_idx >= len(vec):
        return ""
    parts = split_str(vec[dc_idx], ':')
    if level >= len(parts):
        return ""
    return parts[level]

def calculate_osd_requirements(shapes: list[str]) -> int:
    total_osds = 0
    for shape in shapes:
        parts = split_str(shape, ':')
        if not parts:
            continue
        
        # Default to 1 if component unspecified, e.g. "2:2" implies 2 racks, 2 hosts, 1 osd?
        # Actually, let's look at the loop logic in gen_crushmap.
        # It iterates ranges (rack_count, host_count, osd_count).
        # count_spec = get_hierarchy_spec... int(spec)
        
        # Level 0 (DC -> Racks)
        racks = int(parts[0]) if len(parts) > 0 and parts[0] else 0
        # Level 1 (Rack -> Hosts)
        hosts = int(parts[1]) if len(parts) > 1 and parts[1] else 0
        # Level 2 (Host -> OSDs)
        osds = int(parts[2]) if len(parts) > 2 and parts[2] else 0
        
        total_osds += racks * hosts * osds
        
    return total_osds

def validate_osd_count(shapes: list[str], osd_num: int):
    required_osds = calculate_osd_requirements(shapes)
    if osd_num < required_osds:
        raise ValueError(
            f"Insufficient OSDs: Configured osd_num={osd_num} but topology requires {required_osds} "
            f"(shapes={shapes})"
        )



def gen_crushmap():
    root = Node("root", 0)
    
    # Ensure params are lists
    shapes = snakemake.params.dc_shape
    if isinstance(shapes, str):
        shapes = [shapes]
        
    weights = snakemake.params.dc_weight
    if isinstance(weights, str):
        weights = [weights]
        
    # Validate OSD count
    osd_num = int(snakemake.params.osd_num)
    validate_osd_count(shapes, osd_num)

    # Iterate over defined shapes (Datacenters)
    for dc_i in range(len(shapes)):
        dc_str = f"dc{dc_i}"
        dc = Node("dc", 0)
        
        # Level 0: Datacenter -> Rack Count
        count_spec_dc = get_hierarchy_spec(shapes, dc_i, 0)
        rack_count = int(count_spec_dc) if count_spec_dc else 0
        
        weight_spec_dc = get_hierarchy_spec(weights, dc_i, 0)
        dc_weight = resolve_val(weight_spec_dc, dc_i, 0, 4.0)
        
        for rack_i in range(rack_count):
            rack = Node("rack", 0)
            
            # Level 1: Rack -> Host Count
            count_spec_rack = get_hierarchy_spec(shapes, dc_i, 1)
            host_count = int(count_spec_rack) if count_spec_rack else 0
            
            weight_spec_rack = get_hierarchy_spec(weights, dc_i, 1)
            rack_weight = resolve_val(weight_spec_rack, rack_i, dc_i, dc_weight)
            
            for host_i in range(host_count):
                host = Node("host", 0)
                
                # Level 2: Host -> OSD Count
                count_spec_host = get_hierarchy_spec(shapes, dc_i, 2)
                osd_count = int(count_spec_host) if count_spec_host else 0
                
                weight_spec_host = get_hierarchy_spec(weights, dc_i, 2)
                host_weight = resolve_val(weight_spec_host, host_i, rack_i, rack_weight)
                
                for osd_i in range(osd_count):
                    # Resolve OSD weight
                    osd_weight = resolve_val(weight_spec_host, osd_i, host_i, host_weight)
                    
                    osd = Node("osd.", osd_weight)
                    host.add_child(osd)
                
                rack.add_child(host)
            
            dc.add_child(rack)
        
        root.add_child(dc)

    ctx = {"roots": [root]}

    content = template.render(ctx)

    output_file = snakemake.output[0]
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_file, "w") as f:
        f.write(content)

    print(f"Generated crushmap at {output_file}")

gen_crushmap()
