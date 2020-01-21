import sys
import os

import openmc


def build_homog_input_files(order, source_loc, bc, prob_size,
                                  num_neutrons, num_batches, SE_dim,
                                  seed, dir):
    os.system("mkdir -p {dir}".format(dir=dir))
    os.chdir(dir)

    xmin = -1.*prob_size[0]/2.
    xmax = prob_size[0]/2.
    ymin = -1.*prob_size[1]/2.
    ymax = prob_size[1]/2.
    zmin = -1.*prob_size[2]/2.
    zmax = prob_size[2]/2.

    # Instantiate some Materials and register the appropriate Nuclides
    uranyl_sulf = openmc.Material(name='Uranyl Sulfate')
    uranyl_sulf.set_density("atom/b-cm", 9.9035E-02)
    uranyl_sulf.add_nuclide("U235", 9.6795E-05)
    uranyl_sulf.add_nuclide("U234", 7.4257E-07)
    uranyl_sulf.add_nuclide("U238", 5.5518E-04)
    uranyl_sulf.add_nuclide("S32" , 6.5272E-04)
    uranyl_sulf.add_nuclide("O16" , 3.5185E-02)
    uranyl_sulf.add_nuclide("H1"  , 6.2538E-02)
    uranyl_sulf.add_s_alpha_beta('c_H_in_H2O')
    uranyl_sulf.depletable = False

    # Instantiate a Materials collection and export to XML
    materials_file = openmc.Materials([uranyl_sulf])
    materials_file.export_to_xml()

    # Instantiate planar surfaces
    x1 = openmc.XPlane(x0=xmin)
    x2 = openmc.XPlane(x0=xmax)
    y1 = openmc.YPlane(y0=ymin)
    y2 = openmc.YPlane(y0=ymax)
    z1 = openmc.ZPlane(z0=zmin)
    z2 = openmc.ZPlane(z0=zmax)

    # Set boundary conditions
    surface_list = [x1, x2, y1, y2, z1, z2]
    for i in range(len(surface_list)):
        surface = surface_list[i]
        surface.boundary_type = bc[i]

    # Define Cell, Region, and Fill
    uran_sulf_sol = openmc.Cell(name='uranyl sulfate solution')
    uran_sulf_sol.region = +x1 & -x2 & +y1 & -y2 & +z1 & -z2
    uran_sulf_sol.fill = uranyl_sulf

    # Instantiate root universe
    root = openmc.Universe(name='root universe')
    root.add_cells([uran_sulf_sol])

    # Instantiate a Geometry, register the root Universe, and export to XML
    geometry = openmc.Geometry(root)
    geometry.export_to_xml()

    # Define runtime settings
    settings = openmc.Settings()
    point_source = openmc.stats.Point(xyz=source_loc)
    settings.source = openmc.Source(space=point_source)
    settings.batches = num_batches
    settings.inactive = num_batches - 1
    settings.particles = num_neutrons

    # Define entropy mesh
    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left = [xmin, ymin, zmin]
    entropy_mesh.upper_right = [xmax, ymax, zmax]
    entropy_mesh.dimension = SE_dim
    settings.entropy_mesh = entropy_mesh
    settings.seed = seed
    settings.export_to_xml()

    # Create a nu-fission tally
    nu_fiss_tally = openmc.Tally()
    nu_fiss_tally.scores = ['nu-fission']

    # Create a Legendre polynomial expansion filter and add to tally
    expand_filter = openmc.SpatialLegendreFilter(order, 'z', zmin, zmax)
    nu_fiss_tally.filters.append(expand_filter)

    tallies = openmc.Tallies([nu_fiss_tally])
    tallies.export_to_xml()

    os.chdir("./..")


def get_source_locs(prob_type, prob_size = None):
    source_loc_dict = {
        '1d-homog': {
            100.: [(0.0, 0.0, 0.0), (0.0, 0.0, 40.0)],
            200.: [(0.0, 0.0, 0.0), (0.0, 0.0, 90.0)],
            400.: [(0.0, 0.0, 0.0), (0.0, 0.0, 190.0)],
            600.: [(0.0, 0.0, 0.0), (0.0, 0.0, 290.0)],
            800.: [(0.0, 0.0, 0.0), (0.0, 0.0, 390.0)]
        },
        '2d-beavrs':
            [(-161.2773, -161.2773, 220.0, 161.2773, 161.2773, 230.0),
             (0.0, -161.2773, 220.0, 161.2773, 161.2773, 230.0)],
        '3d-exasmr':
            [(-75.26274, -75.26274, 36.007, 75.26274, 75.26274, 236.0066),
             (0.0, -75.26274, 36.007, 75.26274, 75.26274, 236.0066)],
    }

    if prob_size is None:
        return source_loc_dict[prob_type]
    else:
        return source_loc_dict[prob_type][prob_size[2]]


def get_num_batches(prob_type, prob_size = None):
    num_batch_dict = {
        '1d-homog': {
            100.: [399, 499],
            200.: [499, 499],
            400.: [599, 999],
            600.: [999, 2099],
            800.: [1499, 2999]
        },
        '2d-beavrs': [499, 999],
        '3d-exasmr': [149, 199],
    }

    if prob_size is None:
        return num_batch_dict[prob_type]
    else:
        return num_batch_dict[prob_type][prob_size[2]]


def create_file_from_template(replace_dict, template_path, outfile, dir):
    with open(template_path, 'r') as file:
      new_file=file.read()
    for key in replace_dict:
        new_file = new_file.replace(key, str(replace_dict[key]))
    with open(outfile, "w") as file:
      file.write(new_file)
    os.system("mv {} {}".format(outfile, dir))


def build_benchmark_input_files(prob_type, params, dir, num_batch,
                                num_neutrons, seed, source_loc):

    os.system("mkdir -p {dir}".format(dir=dir))

    geometry_replace_dict = {}
    template_path = '../base/geometry_template.xml'
    outfile = 'geometry.xml'
    create_file_from_template(geometry_replace_dict, template_path, outfile,
                              dir)

    materials_replace_dict = {}
    template_path = '../base/materials_template.xml'
    outfile = 'materials.xml'
    create_file_from_template(materials_replace_dict, template_path, outfile,
                              dir)

    SE_dim = params['SE_dim']
    settings_replace_dict = {
        '{num_particles}': str(num_neutrons),
        '{num_batches}': str(num_batch),
        '{num_inactive}': str(num_batch - 1),
        '{source_loc}': ' '.join(str(s) for s in source_loc),
        '{SE_dim}': ' '.join(str(s) for s in SE_dim),
        '{seed}': str(seed)
    }
    template_path = '../base/settings_template.xml'
    outfile = 'settings.xml'
    create_file_from_template(settings_replace_dict, template_path, outfile,
                              dir)

    if prob_type == '2d-beavrs':
        tallies_replace_dict = {
            '{tally_order}': str(params['tally_order'])
        }
        template_path = '../base/tallies_template.xml'
        outfile = 'tallies.xml'
        create_file_from_template(tallies_replace_dict, template_path, outfile,
                                  dir)
    elif prob_type == '3d-exasmr':
        tallies_replace_dict = {
            '{zern_tally_order}': str(params['zern_tally_order']),
            '{leg_tally_order}': str(params['leg_tally_order']),
        }
        template_path = '../base/tallies_template.xml'
        outfile = 'tallies.xml'
        create_file_from_template(tallies_replace_dict, template_path, outfile,
                                  dir)


def generate_input_files(prob_type, seed):
    # Get problem parameters
    params = get_problem_params(prob_type)

    if prob_type == '1d-homog':
        target_dir = '../{}/seed{}'.format(prob_type, seed)
        os.system('mkdir -p {}'.format(target_dir))
        os.chdir(target_dir)

        for num_neutrons in params['num_neutrons']:
            order = params['tally_order']
            bc = params['bc']
            SE_dim = params['SE_dim']
            for prob_size in params['prob_sizes']:
                source_locs = get_source_locs(prob_type,
                                              prob_size=prob_size)
                num_batches = get_num_batches(prob_type,
                                              prob_size=prob_size)

                for i in range(len(source_locs)):
                    source_loc = source_locs[i]
                    num_batch = num_batches[i]
                    dir = "{neutrons}n-{size}cm-{loc}source" \
                        .format(loc=source_loc[2], neutrons=num_neutrons,
                                size=prob_size[2])
                    build_homog_input_files(order, source_loc, bc,
                                            prob_size, num_neutrons,
                                            num_batch, SE_dim, seed,
                                            dir)
    else:
        target_dir = '../{}/seed{}'.format(prob_type, seed)
        os.system('mkdir -p {}'.format(target_dir))
        os.chdir(target_dir)

        for num_neutrons in params['num_neutrons']:
            source_locs = get_source_locs(prob_type)
            num_batches = get_num_batches(prob_type)
            for i in range(len(source_locs)):
                source_loc = source_locs[i]
                num_batch = num_batches[i]
                dir = "{neutrons}n".format(neutrons=num_neutrons)
                if i == 1:
                    dir += "-offset"
                build_benchmark_input_files(prob_type, params, dir, num_batch,
                                            num_neutrons, seed, source_loc)


def get_problem_params(prob_type):
    problem_param_dict = {
        '1d-homog': {
            'tally_order': 30,
            'num_neutrons': [10000, 100000, 1000000, 10000000],
            'SE_dim': [1, 1, 1000],
            'bc': ['reflective', 'reflective', 'reflective', 'reflective',
                   'vacuum', 'vacuum'],
            'prob_sizes': [(10., 10., 100.), (10., 10., 200.),
                           (10., 10., 400.), (10., 10., 600.),
                           (10., 10., 800.)]
        },
        '2d-beavrs': {
            'tally_order': 20,
            'num_neutrons': [100000, 1000000, 4000000, 10000000, 40000000],
            'SE_dim': [68, 68, 1]
        },
        '3d-exasmr': {
            'leg_tally_order': 30,
            'zern_tally_order': 20,
            'num_neutrons': [100000, 1000000, 10000000, 40000000],
            'SE_dim': [28, 28, 20]
        }
    }

    return problem_param_dict[prob_type]


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print('Usage: generate_examples.py [problem_type] [seed #]')
        sys.exit()

    # Get command line arguments
    prob_type = sys.argv[1]
    try:
        seed = int(sys.argv[2])
    except ValueError:
        print('Seed number must be of type int')
        sys.exit()

    # Generate OpenMC and batch script input files
    generate_input_files(prob_type, seed)
