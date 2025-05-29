import re
import openmc
import numpy as np
import os

# Material parsing functions from mcnp_to_openmc.py
def zaid_to_openmc_nuclide(zaid):
    """Convert MCNP ZAID (e.g., 92235.11c) to OpenMC nuclide name (e.g., U235)."""
    zaid = zaid.split('.')[0]
    z = int(zaid[:-3])
    a = int(zaid[-3:])
    element_map = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
        23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
        30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
        37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
        44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
        51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
        58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
        65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
        72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
        79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
        86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
        93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
        100: 'Fm'
    }
    metastable_map = {'95242': 'Am242_m1', '99254': 'Es254_m1'}
    if zaid in metastable_map:
        return metastable_map[zaid]
    element = element_map.get(z, f'Element{z}')
    if element == f'Element{z}':
        print(f"Warning: Element Z={z} not found. Using {element}.")
    return f'{element}{a}'

def parse_mcnp_material_card(filename):
    """Parse MCNP material card and return a list of material definitions."""
    materials = []
    current_material = None
    atom_fractions = {}
    total_density = None
    material_id = None
    with open(filename, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('c'):
                i += 1
                continue
            if re.match(r'^m\d+', line):
                if current_material is not None:
                    materials.append({
                        'id': material_id,
                        'atom_fractions': atom_fractions,
                        'total_density': total_density
                    })
                material_id = int(line.split()[0][1:])
                atom_fractions = {}
                total_density = None
                current_material = line
                i += 1
                while i < len(lines) and not re.match(r'^m\d+', lines[i].strip()) and not lines[i].strip().startswith('c'):
                    line = lines[i].strip()
                    if not line:
                        i += 1
                        continue
                    if 'Total' in line and 'E' in line:
                        try:
                            total_density = float(re.search(r'Total\s+([\d.E+-]+)', line).group(1))
                        except:
                            total_density = None
                    else:
                        tokens = line.split()
                        j = 0
                        while j < len(tokens) - 1:
                            if re.match(r'\d+\.\d+c', tokens[j]):
                                zaid = tokens[j]
                                try:
                                    fraction = float(tokens[j + 1])
                                    atom_fractions[zaid] = atom_fractions.get(zaid, 0) + fraction
                                except (ValueError, IndexError):
                                    print(f"Warning: Invalid fraction for ZAID {zaid} in material m{material_id}")
                                j += 2
                            else:
                                j += 1
                    i += 1
        if current_material is not None:
            materials.append({
                'id': material_id,
                'atom_fractions': atom_fractions,
                'total_density': total_density
            })
    return materials

def generate_openmc_materials(materials):
    """Generate OpenMC material definitions from parsed MCNP materials."""
    openmc_materials = {}
    for mat in materials:
        material = openmc.Material(material_id=mat['id'])
        material.name = f'm{mat["id"]}'
        total_fraction = sum(mat['atom_fractions'].values())
        if total_fraction == 0:
            print(f"Warning: Material m{mat['id']} has no valid nuclides. Skipping.")
            continue
        for zaid, fraction in mat['atom_fractions'].items():
            nuclide = zaid_to_openmc_nuclide(zaid)
            atom_percent = fraction / total_fraction if total_fraction > 0 else fraction
            material.add_nuclide(nuclide, atom_percent, 'ao')
        if mat['total_density'] is not None:
            material.set_density('atom/b-cm', mat['total_density'])
        else:
            material.set_density('atom/b-cm', total_fraction)
        openmc_materials[mat['id']] = material
    return openmc_materials

def write_openmc_material_card(materials, output_filename):
    """Write OpenMC material definitions to a Python file."""
    with open(output_filename, 'w') as f:
        f.write('import openmc\n\n')
        f.write('# OpenMC Material Definitions\n')
        f.write('materials = {}\n\n')
        for mat_id, mat in materials.items():
            f.write(f'# Material m{mat_id}\n')
            f.write(f'materials[{mat_id}] = openmc.Material(material_id={mat_id}, name="m{mat_id}")\n')
            nuclide_densities = mat.get_nuclide_atom_densities()
            for nuclide, value in nuclide_densities.items():
                fraction = value[1] if isinstance(value, tuple) else value
                f.write(f'materials[{mat_id}].add_nuclide("{nuclide}", {fraction}, "ao")\n')
            f.write(f'materials[{mat_id}].set_density("atom/b-cm", {mat.density})\n\n')

def generate_materials():
    """Generate openmc_materials.py from matcard.txt."""
    material_file = 'matcard.txt'
    if not os.path.exists(material_file):
        raise FileNotFoundError(f"Material file '{material_file}' not found")
    try:
        parsed_materials = parse_mcnp_material_card(material_file)
        openmc_materials = generate_openmc_materials(parsed_materials)
        write_openmc_material_card(openmc_materials, 'openmc_materials.py')
        print("Generated openmc_materials.py")
        return openmc_materials
    except Exception as e:
        print(f"Error: Failed to generate materials: {e}")
        return {}

def parse_mcnp_input_file(filename):
    """Parse MCNP input file containing cell and surface cards."""
    surfaces = []
    cells = []
    cell_ids = set()
    surface_ids = set()
    current_cell = None
    fill_array = None
    in_fill_section = False
    line_number = 0
    with open(filename, 'r') as f:
        for line in f:
            line_number += 1
            line = line.strip()
            if not line or line.startswith('c'):
                continue
            tokens = line.split('$')[0].strip().split()
            if not tokens:
                continue
            if re.match(r'^\d+', tokens[0]):
                if current_cell and in_fill_section:
                    current_cell['fill_array'] = fill_array
                    in_fill_section = False
                try:
                    id_num = int(tokens[0])
                    if len(tokens) > 1 and tokens[1].lower() in ['pz', 'cz', 'c/z', 'c/x', 'c/y', 'px', 'py', 'rhp', 'sph', 'so']:
                        if id_num in cell_ids:
                            print(f"Warning: Surface ID {id_num} at line {line_number} conflicts with cell ID")
                            continue
                        if id_num in surface_ids:
                            print(f"Warning: Duplicate surface ID {id_num} at line {line_number}")
                            continue
                        surface_type = tokens[1].lower()
                        params = [float(x) for x in tokens[2:] if x.replace('.', '').replace('-', '').replace('E', '').isdigit()]
                        surfaces.append({
                            'id': id_num,
                            'type': surface_type,
                            'params': params
                        })
                        surface_ids.add(id_num)
                    else:
                        if id_num in surface_ids:
                            print(f"Warning: Cell ID {id_num} at line {line_number} conflicts with surface ID")
                            continue
                        if id_num in cell_ids:
                            print(f"Warning: Duplicate cell ID {id_num} at line {line_number}")
                            continue
                        material_id = int(tokens[1])
                        density = float(tokens[2]) if material_id != 0 else None
                        region_tokens = []
                        universe = None
                        lattice = None
                        fill = None
                        translation = None
                        i = 3
                        while i < len(tokens):
                            if tokens[i].startswith('u='):
                                universe = int(tokens[i].split('=')[1])
                            elif tokens[i].startswith('lat='):
                                lattice = int(tokens[i].split('=')[1])
                            elif tokens[i].startswith('fill'):
                                fill_str = tokens[i].split('=')[1]
                                if '(' in fill_str:
                                    fill_id = int(fill_str.split('(')[0])
                                    trans = fill_str.split('(')[1].strip(')').split()
                                    translation = [float(x) for x in trans]
                                    fill = (fill_id, translation)
                                else:
                                    fill = int(fill_str.split(':')[0])
                            elif tokens[i].startswith('IMP:n='):
                                pass  # Ignored
                            else:
                                region_tokens.append(tokens[i])
                            i += 1
                        region_str = ' '.join(region_tokens)
                        current_cell = {
                            'id': id_num,
                            'material_id': material_id,
                            'density': density,
                            'region': region_str,
                            'universe': universe,
                            'lattice': lattice,
                            'fill': fill
                        }
                        cells.append(current_cell)
                        cell_ids.add(id_num)
                except (ValueError, IndexError) as e:
                    print(f"Warning: Failed to parse line {line_number}: '{line}': {e}")
            elif current_cell and re.match(r'^\d+\s+\d+', line):
                in_fill_section = True
                if fill_array is None:
                    fill_array = []
                try:
                    fill_array.append([int(x) for x in tokens if x.isdigit()])
                except ValueError:
                    print(f"Warning: Invalid fill array at line {line_number}: '{line}'")
    if current_cell and in_fill_section:
        current_cell['fill_array'] = fill_array
    universes = {}
    lattices = {}
    for cell in cells:
        u = cell['universe'] if cell['universe'] is not None else 0
        universes.setdefault(u, []).append(cell)
    for cell in cells:
        if cell['lattice']:
            lattices[cell['id']] = cell
    return surfaces, cells, universes, lattices

def generate_openmc_surfaces(surfaces):
    """Generate OpenMC surface definitions."""
    openmc_surfaces = {}
    has_hexagonal_prism = hasattr(openmc, 'HexagonalPrism')
    for surf in surfaces:
        try:
            surface_id = surf['id']
            surface_type = surf['type']
            params = surf['params']
            if surface_type == 'pz':
                surface = openmc.ZPlane(z0=params[0], name=f's{surface_id}')
            elif surface_type == 'cz':
                surface = openmc.ZCylinder(r=params[0], name=f's{surface_id}')
            elif surface_type == 'c/z':
                surface = openmc.ZCylinder(x0=params[0], y0=params[1], r=params[2], name=f's{surface_id}')
            elif surface_type == 'c/x':
                surface = openmc.XCylinder(y0=params[0], z0=params[1], r=params[2], name=f's{surface_id}')
            elif surface_type == 'c/y':
                surface = openmc.YCylinder(x0=params[0], z0=params[1], r=params[2], name=f's{surface_id}')
            elif surface_type == 'px':
                surface = openmc.XPlane(x0=params[0], name=f's{surface_id}')
            elif surface_type == 'py':
                surface = openmc.YPlane(y0=params[0], name=f's{surface_id}')
            elif surface_type == 'rhp':
                if has_hexagonal_prism:
                    x0, y0, z0, dx, dy, dz, edge_length = params[:7]
                    center = (x0 + dx/2, y0 + dy/2)
                    surface = openmc.HexagonalPrism(center=center, edge_length=edge_length, orientation='y', name=f's{surface_id}')
                else:
                    print(f"Warning: HexagonalPrism not available in OpenMC. Approximating rhp surface {surface_id} with planes.")
                    x0, y0, z0, dx, dy, dz, edge_length = params[:7]
                    center_x, center_y = x0 + dx/2, y0 + dy/2
                    # Define six planes for a hexagonal prism (y-orientation)
                    r = edge_length  # Distance from center to vertex
                    angles = np.linspace(0, 2 * np.pi, 7)[:-1]  # 6 angles (0, 60, 120, ...)
                    planes = []
                    for i, theta in enumerate(angles):
                        # Plane normal: (cos(theta), sin(theta), 0)
                        # Distance: r * cos(angle between normal and point)
                        nx = np.cos(theta)
                        ny = np.sin(theta)
                        # Point on plane: center + r * (cos(theta + 90), sin(theta + 90))
                        px = center_x + r * np.cos(theta + np.pi/2)
                        py = center_y + r * np.sin(theta + np.pi/2)
                        # Plane equation: nx*x + ny*y = nx*px + ny*py
                        d = nx * px + ny * py
                        plane = openmc.Plane(a=nx, b=ny, c=0, d=d, name=f's{surface_id}_p{i}')
                        planes.append(plane)
                    # Combine planes: inside all planes
                    surface = planes[0]
                    for p in planes[1:]:
                        surface &= p
            elif surface_type == 'sph':
                surface = openmc.Sphere(x0=params[0], y0=params[1], z0=params[2], r=params[3], name=f's{surface_id}')
            elif surface_type == 'so':
                surface = openmc.Sphere(r=params[0], name=f's{surface_id}')
            else:
                print(f"Warning: Unsupported surface type '{surface_type}' for surface {surface_id}")
                continue
            openmc_surfaces[surface_id] = surface
        except (IndexError, ValueError) as e:
            print(f"Warning: Failed to create surface {surface_id}: {e}")
    return openmc_surfaces

def parse_mcnp_region(region_str, surfaces):
    """Parse MCNP region specification into OpenMC region."""
    if not region_str:
        return None
    region_str = region_str.replace('(', '(').replace(')', ')').replace(':', '|')
    tokens = region_str.split()
    stack = []
    current = []
    for token in tokens:
        try:
            if token.isdigit() or (token.startswith('-') and token[1:].isdigit()):
                surf_id = int(token)
                if surf_id in surfaces:
                    current.append(-surfaces[surf_id] if token.startswith('-') else +surfaces[surf_id])
                else:
                    print(f"Warning: Surface {surf_id} not found in region '{region_str}'")
            elif token == '(':
                stack.append(current)
                current = []
            elif token == ')':
                if current:
                    region = current[0]
                    for r in current[1:]:
                        region &= r
                    current = stack.pop()
                    current.append(region)
            elif token == '|':
                if len(current) >= 2:
                    region = current[-2] | current[-1]
                    current = current[:-2] + [region]
        except Exception as e:
            print(f"Warning: Error parsing region token '{token}' in '{region_str}': {e}")
    if current:
        region = current[0]
        for r in current[1:]:
            region &= r
        return region
    return None

def generate_openmc_geometry(cells, universes, lattices, surfaces, materials):
    """Generate OpenMC cells, universes, and lattices."""
    openmc_cells = {}
    openmc_universes = {}
    openmc_lattices = {}
    for cell in cells:
        try:
            cell_id = cell['id']
            region = parse_mcnp_region(cell['region'], surfaces)
            material_id = cell['material_id']
            # Handle m278 as m78 (possible typo)
            if material_id == 278:
                print(f"Warning: Material m278 in cell {cell_id} not found in matcard.txt. Assuming m78 (YH).")
                material_id = 78
            if material_id == 0:
                material = None
            else:
                if material_id in materials:
                    material = materials[material_id]
                    if cell['density'] is not None:
                        density = abs(cell['density']) if cell['density'] < 0 else cell['density']
                        unit = 'g/cm^3'
                        if cell['density'] < 0:
                            print(f"Info: Cell {cell_id} uses MCNP density {density} g/cm³, overriding matcard.txt")
                        elif cell['density'] > 0:
                            print(f"Warning: Cell {cell_id} has positive density {density}. Assuming g/cm³.")
                        material.set_density(unit, density)
                    else:
                        print(f"Warning: No density specified for cell {cell_id}. Using matcard.txt atom/b-cm.")
                else:
                    print(f"Error: Material m{material_id} not found for cell {cell_id}. Setting to void.")
                    material = None
            fill = None
            translation = None
            if cell['fill']:
                if isinstance(cell['fill'], tuple):
                    fill_id, translation = cell['fill']
                    fill = openmc_universes.get(fill_id)
                    if fill is None:
                        print(f"Warning: Universe {fill_id} not found for fill in cell {cell_id}")
                else:
                    fill = openmc_universes.get(cell['fill'])
                    if fill is None:
                        print(f"Warning: Universe {cell['fill']} not found for fill in cell {cell_id}")
            openmc_cell = openmc.Cell(cell_id=cell_id, name=f'c{cell_id}', region=region, fill=fill, translation=translation)
            if material:
                openmc_cell.fill = material
            openmc_cells[cell_id] = openmc_cell
        except Exception as e:
            print(f"Warning: Failed to create cell {cell_id}: {e}")
    for u_id, u_cells in universes.items():
        try:
            universe_cells = [openmc_cells[c['id']] for c in u_cells if c['id'] in openmc_cells]
            openmc_universes[u_id] = openmc.Universe(cells=universe_cells, name=f'u{u_id}')
        except Exception as e:
            print(f"Warning: Failed to create universe {u_id}: {e}")
    for lat_id, lat_cell in lattices.items():
        try:
            if lat_cell['lattice'] == 2:
                fill_array = lat_cell['fill_array']
                lattice = openmc.HexLattice(name=f'lat{lat_id}')
                lattice.center = (0, 0)
                lattice.pitch = (2.54524866 * 2 / np.sqrt(3),) if not fill_array else (len(fill_array[0]) * 2 / np.sqrt(3),)
                lattice.orientation = 'y'
                universe_ids = []
                for row in fill_array:
                    row_universes = []
                    for u_id in row:
                        row_universes.append(openmc_universes.get(u_id, None))
                        if u_id not in openmc_universes:
                            print(f"Warning: Universe {u_id} not found for lattice {lat_id}")
                    universe_ids.append(row_universes)
                lattice.universes = universe_ids[::-1]
                openmc_lattices[lat_id] = lattice
                openmc_cells[lat_id].fill = lattice
            elif lat_cell['lattice'] == 1:
                print(f"Warning: Rectangular lattice (lat=1) for cell {lat_id} not fully supported")
            else:
                print(f"Warning: Unsupported lattice type {lat_cell['lattice']} for cell {lat_id}")
        except Exception as e:
            print(f"Warning: Failed to create lattice for cell {lat_id}: {e}")
    for cell in cells:
        try:
            if cell['id'] in openmc_cells:
                u_id = cell['universe'] if cell['universe'] is not None else 0
                if u_id in openmc_universes:
                    openmc_cells[cell['id']].universe = openmc_universes[u_id]
                else:
                    print(f"Warning: Universe {u_id} not found for cell {cell['id']}")
        except Exception as e:
            print(f"Warning: Failed to set universe for cell {cell['id']}: {e}")
    return openmc_cells, openmc_universes, openmc_lattices

def write_openmc_geometry(cells, universes, lattices, surfaces, output_filename):
    """Write OpenMC geometry definitions to a Python file."""
    try:
        with open(output_filename, 'w') as f:
            f.write('import openmc\n')
            f.write('from openmc_materials import materials\n\n')
            f.write('# OpenMC Geometry Definitions\n')
            f.write('# Surfaces\n')
            for surf_id, surf in surfaces.items():
                f.write(f's{surf_id} = {repr(surf)}\n')
            f.write('\n# Cells\n')
            for cell_id, cell in cells.items():
                f.write(f'c{cell_id} = {repr(cell)}\n')
            f.write('\n# Universes\n')
            for u_id, universe in universes.items():
                f.write(f'u{u_id} = {repr(universe)}\n')
            f.write('\n# Lattices\n')
            for lat_id, lattice in lattices.items():
                f.write(f'lat{lat_id} = {repr(lattice)}\n')
            f.write('\n# Geometry\n')
            f.write('geometry = openmc.Geometry()\n')
            f.write('geometry.root_universe = u0\n')
            f.write('geometry.export_to_xml()\n')
    except Exception as e:
        print(f"Error: Failed to write geometry file: {e}")

def main():
    geometry_file = 'mcnp_input.txt'
    if not os.path.exists(geometry_file):
        raise FileNotFoundError(f"Geometry file '{geometry_file}' not found")
    materials = generate_materials()
    try:
        surfaces, cells, universes, lattices = parse_mcnp_input_file(geometry_file)
    except Exception as e:
        print(f"Error: Failed to parse '{geometry_file}': {e}")
        return
    try:
        openmc_surfaces = generate_openmc_surfaces(surfaces)
        openmc_cells, openmc_universes, openmc_lattices = generate_openmc_geometry(
            cells, universes, lattices, openmc_surfaces, materials
        )
    except Exception as e:
        print(f"Error: Failed to generate OpenMC geometry: {e}")
        return
    try:
        write_openmc_geometry(openmc_cells, openmc_universes, openmc_lattices, openmc_surfaces, 'openmc_geometry.py')
        print("Generated openmc_geometry.py")
    except Exception as e:
        print(f"Error: Failed to write output: {e}")

if __name__ == '__main__':
    main()