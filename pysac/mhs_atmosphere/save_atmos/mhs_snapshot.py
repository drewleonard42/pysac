# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 17:56:58 2014

@author: sm1fg
"""
import numpy as np
import pysac.io.gdf_writer as gdf
import h5py
import astropy.units as u

##============================================================================
## Save a file!!!
##============================================================================
""" For large data arrays this has been producing overfull memory so moved to
dedicated serial function, which should be parallel and moved ahead of gather
for plotting -- maybe handle plotting from hdf5 also

This file is potentially large - recommended to mkdir hdf5 in /data/${USER}
and add symlink to ${HOME} to avoid exceeding quota.
"""
def save_SACvariables(
                    filename,
                    rho,
                    Bx,
                    By,
                    Bz,
                    energy,
                    option_pars,
                    physical_constants,
                    coords,
                    Nxyz,
                    xindex,
                    rank=0,
                    collective=False,
                     ):

    """ Save the background variables for a SAC model in hdf5 (gdf default)
    format after collating the data from mpi sub processes if necessary.
    """
    if rank == 0:
        print'writing',filename
        print'SAC background atmosphere'
    grid_dimensions = [Nxyz[0], Nxyz[1], Nxyz[2]]
    left_edge =  u.Quantity([coords['xmin'],
                             coords['ymin'],
                             coords['zmin']]).to(u.m)
    right_edge = u.Quantity([coords['xmax'],
                             coords['ymax'],
                             coords['zmax']]).to(u.m)

    gamma = 5./3.
    dummy = np.zeros(rho.shape)
    simulation_parameters = gdf.SimulationParameters([
                        ['boundary_conditions', np.zeros(6) + 2],
                        ['cosmological_simulation', 0          ],
                        ['current_iteration', 0                ],
                        ['current_time', 0.0                   ],
                        ['dimensionality', 3                   ],
                        ['domain_dimensions', grid_dimensions  ],
                        ['domain_left_edge', left_edge         ],
                        ['domain_right_edge', right_edge       ],
                        ['eta', 0.0                            ],
                        ['field_ordering', 0                   ],
                        ['gamma', gamma                   ],
                        ['gravity0', 0.0                       ],
                        ['gravity1', 0.0                       ],
                        ['gravity2', physical_constants['gravity'] ],
                        ['nu', 0.0                             ],
                        ['num_ghost_zones', 0                  ],
                        ['refine_by', 0                        ],
                        ['unique_identifier', 'sacgdf2014'     ]
                        ])
    if collective:
        from mpi4py import MPI
        gdf_file = gdf.create_file(h5py.File(filename,'w'),
                                   simulation_parameters, grid_dimensions,
                                   driver='mpio', comm=MPI.COMM_WORLD)
    else:
        gdf_file = gdf.create_file(h5py.File(filename,'w'), simulation_parameters, grid_dimensions)
    
    gdf.write_field(gdf_file, rho,
                          'density_bg',
                          'Background Density'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('kg/m**3'),
                          'density_pert',
                          'Perturbation Density'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          energy,
                          'internal_energy_bg',
                          'Background Internal Energy'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('kg/(m s**2)'),
                          'internal_energy_pert',
                          'Perturbation Internal Energy'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, Bx,
                          'mag_field_x_bg',
                          'x Component of Background Magnetic Field'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('T'),
                          'mag_field_x_pert',
                          'x Component of Pertubation Magnetic Field'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, By,
                          'mag_field_y_bg',
                          'y Component of Background Magnetic Field'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('T'),
                          'mag_field_y_pert',
                          'y Component of Pertubation Magnetic Field'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, Bz,
                          'mag_field_z_bg',
                          'z Component of Background Magnetic Field'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('T'),
                          'mag_field_z_pert',
                          'z Component of Pertubation Magnetic Field'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('m/s'),
                          'velocity_x',
                          'x Component of Velocity'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('m/s'),
                          'velocity_y',
                          'y Component of Velocity'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file, dummy*u.Unit('m/s'),
                          'velocity_z',
                          'z Component of Velocity'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )

    gdf_file.close()
    if collective:
        mpi_save_SACvariables(
                              filename,
                              ((rho,   'density_bg'        ),
                               (energy,'internal_energy_bg'),
                               (Bx,    'mag_field_x_bg'    ),
                               (By,    'mag_field_y_bg'    ),
                               (Bz,    'mag_field_z_bg'    )),
                              xindex,
                             )


#============================================================================

def save_SACsources(
                    sourcesfile,
                    Fx,
                    Fy,
                    option_pars,
                    physical_constants,
                    coords,
                    Nxyz,
                    xindex,
                    rank=0,
                    collective=False,
                   ):
    """ Save the balancing forces for a SAC model with multiple flux tubes in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    if rank == 0:
        print'writing',sourcesfile
        print'SAC background source terms'

    grid_dimensions = [Nxyz[0], Nxyz[1], Nxyz[2]]
    left_edge =  u.Quantity([coords['xmin'],
                             coords['ymin'],
                             coords['zmin']]).to(u.m)
    right_edge = u.Quantity([coords['xmax'],
                             coords['ymax'],
                             coords['zmax']]).to(u.m)

    simulation_parameters = gdf.SimulationParameters([
                        ['boundary_conditions', np.zeros(6) + 2],
                        ['cosmological_simulation', 0          ],
                        ['current_iteration', 0                ],
                        ['current_time', 0.0                   ],
                        ['dimensionality', 3                   ],
                        ['domain_dimensions', grid_dimensions  ],
                        ['domain_left_edge', left_edge         ],
                        ['domain_right_edge', right_edge       ],
                        ['eta', 0.0                            ],
                        ['field_ordering', 0                   ],
                        ['gamma', physical_constants['gamma']  ],
                        ['gravity0', 0.0                       ],
                        ['gravity1', 0.0                       ],
                        ['gravity2', physical_constants['gravity']],
                        ['nu', 0.0                             ],
                        ['num_ghost_zones', 0                  ],
                        ['refine_by', 0                        ],
                        ['unique_identifier', 'sacgdf2014'     ]
                        ])
    if collective:
        from mpi4py import MPI
        gdf_file = gdf.create_file(h5py.File(sourcesfile,'w'),
                                   simulation_parameters, grid_dimensions,
                                   driver='mpio', comm=MPI.COMM_WORLD)
    else:
        gdf_file = gdf.create_file(h5py.File(sourcesfile,'w'),
                                   simulation_parameters, grid_dimensions)

    gdf.write_field(gdf_file,
                          Fx,
                          'balancing_force_x_bg',
                          'x Component of Background Balancing Force'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          Fy,
                          'balancing_force_y_bg',
                          'y Component of Background Balancing Force'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf_file.close()
    if collective:
        mpi_save_SACvariables(
                              sourcesfile,
                              ((Fx,'balancing_force_x_bg'),
                               (Fy,'balancing_force_y_bg')),
                              xindex,
                              )

#============================================================================

def save_auxilliary3D(
                    auxfile,
                    pressure_m,
                    rho_m,
                    temperature,
                    pbeta,
                    alfven,
                    cspeed,
                    dxB2,
                    dyB2,
                    option_pars,
                    physical_constants,
                    coords,
                    Nxyz,
                    xindex,
                    rank=0,
                    collective=False,
                   ):
    """ Save auxilliary variables for use in plotting background setup in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    if rank == 0:
        print'writing',auxfile
        print'non-SAC 3D auxilliary data for plotting'

    grid_dimensions = [Nxyz[0], Nxyz[1], Nxyz[2]]
    left_edge =  u.Quantity([coords['xmin'],
                             coords['ymin'],
                             coords['zmin']]).to(u.m)
    right_edge = u.Quantity([coords['xmax'],
                             coords['ymax'],
                             coords['zmax']]).to(u.m)
    

    simulation_parameters = gdf.SimulationParameters([
                        ['boundary_conditions', np.zeros(6) + 2],
                        ['cosmological_simulation', 0          ],
                        ['current_iteration', 0                ],
                        ['current_time', 0.0                   ],
                        ['dimensionality', 3                   ],
                        ['domain_dimensions', grid_dimensions  ],
                        ['domain_left_edge', left_edge         ],
                        ['domain_right_edge', right_edge       ],
                        ['eta', 0.0                            ],
                        ['field_ordering', 0                   ],
                        ['gamma', physical_constants['gamma']  ],
                        ['gravity0', 0.0                       ],
                        ['gravity1', 0.0                       ],
                        ['gravity2', physical_constants['gravity']],
                        ['nu', 0.0                             ],
                        ['num_ghost_zones', 0                  ],
                        ['refine_by', 0                        ],
                        ['unique_identifier', 'sacgdf2014'     ]
                        ])
    if collective:
        from mpi4py import MPI
        gdf_file = gdf.create_file(h5py.File(auxfile,'w'),
                                   simulation_parameters, grid_dimensions,
                                   driver='mpio', comm=MPI.COMM_WORLD)
    else:
        gdf_file = gdf.create_file(h5py.File(auxfile,'w'),
                                   simulation_parameters, grid_dimensions)

    gdf.write_field(gdf_file,
                          pressure_m,
                          'pressure_mhs',
                          'Background magneto-pressure balance'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          rho_m,
                          'density_mhs',
                          'Background magneto-density balance'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          temperature,
                          'temperature',
                          'Background temperature'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          pbeta,
                          'plasma_beta',
                          'Background plasma beta'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          alfven,
                          'alfven_speed',
                          'Background Alfven speed'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          cspeed,
                          'sound_speed',
                          'Background sound speed'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          dxB2,
                          'mag_tension_x',
                          'x-component background magnetic tension'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          dyB2,
                          'mag_tension_y',
                          'y-component background magnetic tension'
                          ,field_shape=grid_dimensions
                          ,arr_slice=xindex
                          ,collective=collective
                          )
    gdf_file.close()
    if collective:
        mpi_save_SACvariables(
                              auxfile,
                              ((pressure_m,
                                'pressure_mhs'),
                               ( rho_m,
                                'density_mhs'),
                               (temperature,
                                'temperature'),
                               (pbeta,
                                'plasma_beta'),
                               (alfven,
                                'alfven_speed'),
                               (cspeed,
                                'sound_speed'),
                               (dxB2,
                                'mag_tension_x'),
                               (dyB2,
                                'mag_tension_y')),
                              xindex,
                              )

#============================================================================

def save_auxilliary1D(
                    auxfile,
                    pressure_Z,
                    rho_Z,
                    Rgas_Z,
                    option_pars,
                    physical_constants,
                    coords,
                    Nxyz,
                    rank=0,
                    collective=False,
                   ):
    """ Save auxilliary variables for use in plotting background setup in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    if rank == 0:
        print'writing',auxfile
        print'non-SAC 1D auxilliary data for plotting'

    grid_dimensions = [2, 2, Nxyz[2]] #dims > 1 to be read by yt
    left_edge =  u.Quantity([coords['xmin'],
                             coords['ymin'],
                             coords['zmin']]).to(u.m)
    right_edge = u.Quantity([coords['xmax'],
                             coords['ymax'],
                             coords['zmax']]).to(u.m)
    
    pressureHS = u.Quantity(np.zeros(grid_dimensions),
                            unit=pressure_Z.to("Pa").unit)
    rhoHS = u.Quantity(np.zeros(grid_dimensions),
                       unit=rho_Z.to('kg/m**3').unit)
    RgasHS = u.Quantity(np.zeros(grid_dimensions),
                        unit=Rgas_Z.to('m**2/(K s**2)').unit)
    pressureHS[:] = pressure_Z
    rhoHS[:] = rho_Z
    RgasHS[:] = Rgas_Z
    simulation_parameters = gdf.SimulationParameters([
                            ['boundary_conditions', np.zeros(6) + 2],
                            ['cosmological_simulation', 0          ],
                            ['current_iteration', 0                ],
                            ['current_time', 0.0                   ],
                            ['dimensionality', 3                   ],
                            ['domain_dimensions', grid_dimensions  ],
                            ['domain_left_edge', left_edge         ],
                            ['domain_right_edge', right_edge       ],
                            ['eta', 0.0                            ],
                            ['field_ordering', 0                   ],
                            ['gamma', physical_constants['gamma']  ],
                            ['gravity0', 0.0                       ],
                            ['gravity1', 0.0                       ],
                            ['gravity2', physical_constants['gravity']],
                            ['nu', 0.0                             ],
                            ['num_ghost_zones', 0                  ],
                            ['refine_by', 0                        ],
                            ['unique_identifier', 'sacgdf2014'     ]
                            ])

#    import pdb; pdb.set_trace()
    gdf_file = gdf.create_file(h5py.File(auxfile,'w'), simulation_parameters, grid_dimensions)

    gdf.write_field(gdf_file,
                          pressureHS,
                          'pressure_HS',
                          'Background 1D hydrostatic-pressure'
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          rhoHS,
                          'density_HS',
                          'Background 1D hydrostatic-density'
                          ,collective=collective
                          )
    gdf.write_field(gdf_file,
                          RgasHS,
                          'ideal_gas_constant_HS',
                          'Background 1D hydrostatic-R_gas'
                          ,collective=collective
                          )
    gdf_file.close()
    

#=============================================================================
## Save a file!!!
##============================================================================
""" For large data arrays this has been producing overfull memory so moved to
dedicated serial function, which should be parallel and moved ahead of gather
for plotting -- maybe handle plotting from hdf5 also

This file is potentially large - recommended to mkdir hdf5 in /data/${USER}
and add symlink to ${HOME} to avoid exceeding quota.
"""
def mpi_save_SACvariables(
                    filename,
                    fields,
                    xindex,
                     ):

    """ Save the background variables for a SAC model in hdf5 (gdf default)
    format after collating the data from mpi sub processes if necessary.
    """
    from mpi4py import MPI

    comm = MPI.COMM_WORLD

    ds = h5py.File(filename, 'a', driver='mpio', comm=comm)
    for var, field in fields:
        ds['data']['grid_0000000000'][field][xindex,...] = var
    ds.close()
