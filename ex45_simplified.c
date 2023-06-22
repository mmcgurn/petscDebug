static char help[] = "Simplified example for conduction";

#include <petscdmplex.h>
#include <petscds.h>
#include <petscts.h>

/*
The heat transfer equation assuming constant density, specific heat and, conductivity.  There are no source/sink terms.

    \[\rho C \frac{\partial u}{\partial dt} - k \Delta u = 0\]

In weak form:

   \[ \int_{\Omega} \left[w \rho C \frac{\partial T}{\partial t} + k \frac{\partial w}{\partial x_i} \frac{\partial T}{\partial x_j}  \right] d\Omega + \int_{\Gamma} q w d\Gamma \]

Therefore the integrand for the test function w term is

   \[\rho C \frac{\partial T}{\partial t}\]

The integrand for the test function gradient term w is

   \[k \frac{\partial T}{\partial x_j} \]

The gradient \(\frac{\partial F}{\partial T} \)
for g3 - integrand for the test function gradient and basis function gradient term is \(k \) for the diagonal terms. All others are zero.

The gradient \(\frac{\partial F}{\partial \partial T/\partial t} \) for g0 - integrand for the test and basis function term is \(\rho C \).  This uses the u\_tShift scalar in petsc.  All others are zero.
*/


/*
 * TODO:
 * create a box mesh so I can the dimensions
 * move the flags into hard code
 * figure out how to change the boundary condition
 */
typedef enum {
    specificHeat, conductivity, density, total
} ConductionProperties;

/**
 * Compute the test function integrated.  Note there is only a single field.
 */
static void
wIntegrandTestFunction(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                       const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                       const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants,
                       const PetscScalar constants[], PetscScalar f0[]) {
    f0[0] = constants[density] * constants[specificHeat] * u_t[0];
}

/**
 * Compute the test function integrated.  Note there is only a single field.
 */
static void wIntegrandTestGradientFunction(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[],
                                           const PetscInt uOff_x[],
                                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                           const PetscInt aOff[],
                                           const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                                           const PetscScalar a_x[], PetscReal t, const PetscReal x[],
                                           PetscInt numConstants,
                                           const PetscScalar constants[], PetscScalar f1[]) {
    for (PetscInt d = 0; d < dim; ++d) {
        f1[d] = constants[conductivity] * u_x[d];
    }
}

/**
 * Compute the jacobian term g0 - integrand for the test and basis function term i
 */
static void jacobianG0Term(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[],
                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[],
                           const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                           const PetscScalar a_x[],
                           PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
                           const PetscScalar constants[], PetscScalar g3[]) {
    for (PetscInt d = 0; d < dim; ++d) { g3[d * dim + d] = u_tShift * constants[density] * constants[specificHeat]; }
}


/**
 * Compute the jacobian term g3 - integrand for the test function gradient and basis function gradient term
 */
static void jacobianG3Term(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[],
                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[],
                           const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                           const PetscScalar a_x[],
                           PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
                           const PetscScalar constants[], PetscScalar g3[]) {
    for (PetscInt d = 0; d < dim; ++d) { g3[d * dim + d] = constants[conductivity]; }
}

static PetscErrorCode
essential_boundary_condition(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u,
                             void *ctx) {
    *u = 500.0 + x[0];
    return PETSC_SUCCESS;
}


static PetscErrorCode CreateMesh(MPI_Comm comm, DM *dm) {
    PetscFunctionBeginUser;
    PetscCall(DMCreate(comm, dm));
    PetscCall(DMSetType(*dm, DMPLEX));
    PetscCall(DMSetFromOptions(*dm));
    PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetupProblem(DM dm) {
    PetscDS ds;
    DMLabel label;
    const PetscInt id = 1;

    PetscFunctionBeginUser;
    PetscCall(DMGetLabel(dm, "marker", &label));
    PetscCall(DMGetDS(dm, &ds));
    PetscCall(PetscDSSetJacobian(ds, 0, 0, jacobianG0Term, NULL, NULL, jacobianG3Term));
    PetscCall(PetscDSSetResidual(ds, 0, wIntegrandTestFunction, wIntegrandTestGradientFunction));

    // hard code a boundary condition for now
    PetscCall(DMAddBoundary(dm, DM_BC_ESSENTIAL, "wall", label, 1, &id, 0, 0, NULL,
                            (void (*)(void)) essential_boundary_condition, NULL,
                            NULL, NULL));


    // Set the constant values for this problem
    PetscReal parameterArray[total] = {1000, 1, 1};
    PetscCall(PetscDSSetConstants(ds, total, parameterArray));

    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetupDiscretization(DM dm) {
    DM cdm = dm;
    PetscFE fe;
    DMPolytopeType ct;
    PetscInt dim, cStart;

    PetscFunctionBeginUser;
    PetscCall(DMGetDimension(dm, &dim));
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, NULL));
    PetscCall(DMPlexGetCellType(dm, cStart, &ct));
    /* Create finite element */
    PetscCall(PetscFECreateByCell(PETSC_COMM_SELF, dim, 1, ct, "temp_", -1, &fe));
    PetscCall(PetscObjectSetName((PetscObject) fe, "temperature"));
    /* Set discretization and boundary conditions for each mesh */
    PetscCall(DMSetField(dm, 0, NULL, (PetscObject) fe));
    PetscCall(DMCreateDS(dm));

    // Setup the problem
    PetscCall(SetupProblem(dm));

    // copy over the discratation
    while (cdm) {
        PetscCall(DMCopyDisc(dm, cdm));
        PetscCall(DMGetCoarseDM(cdm, &cdm));
    }
    PetscCall(PetscFEDestroy(&fe));
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetInitialConditions(TS ts, Vec u) {
    DM dm;
    PetscReal t;

    PetscFunctionBeginUser;
    PetscCall(TSGetDM(ts, &dm));
    PetscCall(TSGetTime(ts, &t));

    // Set the initial condition
    VecSet(u, 300.0);

    PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
    DM dm;
    TS ts;
    Vec u;

    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, NULL, help));
    PetscCall(CreateMesh(PETSC_COMM_WORLD, &dm));
    PetscCall(SetupDiscretization(dm));

    PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
    PetscCall(TSSetDM(ts, dm));
    PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, NULL));
    PetscCall(DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, NULL));
    PetscCall(DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, NULL));
    PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSSetFromOptions(ts));
    PetscCall(TSSetComputeInitialCondition(ts, SetInitialConditions));

    PetscCall(DMCreateGlobalVector(dm, &u));
    PetscCall(DMTSCheckFromOptions(ts, u));
    PetscCall(SetInitialConditions(ts, u));
    PetscCall(PetscObjectSetName((PetscObject) u, "temperature"));
    PetscCall(TSSolve(ts, u));

    PetscCall(VecDestroy(&u));
    PetscCall(TSDestroy(&ts));
    PetscCall(DMDestroy(&dm));
    PetscCall(PetscFinalize());
    return 0;
}

/*TEST

  test:
    suffix: 2d_p1
    requires: triangle
    args: -sol_type quadratic_linear -dm_refine 1 -temp_petscspace_degree 1 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    suffix: 2d_p1_exp
    requires: triangle
    args: -sol_type quadratic_linear -dm_refine 1 -temp_petscspace_degree 1 -explicit \
          -ts_type euler -ts_max_steps 4 -ts_dt 1e-3 -ts_monitor_error
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [1.9]
    suffix: 2d_p1_sconv
    requires: triangle
    args: -sol_type quadratic_linear -temp_petscspace_degree 1 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [2.1]
    suffix: 2d_p1_sconv_2
    requires: triangle
    args: -sol_type trig_trig -temp_petscspace_degree 1 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 1e-6 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 4 -convest_num_refine 3 get L_2 convergence rate: [1.2]
    suffix: 2d_p1_tconv
    requires: triangle
    args: -sol_type quadratic_trig -temp_petscspace_degree 1 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 6 -convest_num_refine 3 get L_2 convergence rate: [1.0]
    suffix: 2d_p1_tconv_2
    requires: triangle
    args: -sol_type trig_trig -temp_petscspace_degree 1 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # The L_2 convergence rate cannot be seen since stability of the explicit integrator requires that is be more accurate than the grid
    suffix: 2d_p1_exp_tconv_2
    requires: triangle
    args: -sol_type trig_trig -temp_petscspace_degree 1 -explicit -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type euler -ts_max_steps 4 -ts_dt 1e-4 -lumped 0 -mass_pc_type lu
  test:
    # The L_2 convergence rate cannot be seen since stability of the explicit integrator requires that is be more accurate than the grid
    suffix: 2d_p1_exp_tconv_2_lump
    requires: triangle
    args: -sol_type trig_trig -temp_petscspace_degree 1 -explicit -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type euler -ts_max_steps 4 -ts_dt 1e-4
  test:
    suffix: 2d_p2
    requires: triangle
    args: -sol_type quadratic_linear -dm_refine 0 -temp_petscspace_degree 2 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [2.9]
    suffix: 2d_p2_sconv
    requires: triangle
    args: -sol_type trig_linear -temp_petscspace_degree 2 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00000001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [3.1]
    suffix: 2d_p2_sconv_2
    requires: triangle
    args: -sol_type trig_trig -temp_petscspace_degree 2 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00000001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 3 -convest_num_refine 3 get L_2 convergence rate: [1.0]
    suffix: 2d_p2_tconv
    requires: triangle
    args: -sol_type quadratic_trig -temp_petscspace_degree 2 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 3 -convest_num_refine 3 get L_2 convergence rate: [1.0]
    suffix: 2d_p2_tconv_2
    requires: triangle
    args: -sol_type trig_trig -temp_petscspace_degree 2 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    suffix: 2d_q1
    args: -sol_type quadratic_linear -dm_plex_simplex 0 -dm_refine 1 -temp_petscspace_degree 1 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [1.9]
    suffix: 2d_q1_sconv
    args: -sol_type quadratic_linear -dm_plex_simplex 0 -temp_petscspace_degree 1 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 4 -convest_num_refine 3 get L_2 convergence rate: [1.2]
    suffix: 2d_q1_tconv
    args: -sol_type quadratic_trig -dm_plex_simplex 0 -temp_petscspace_degree 1 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    suffix: 2d_q2
    args: -sol_type quadratic_linear -dm_plex_simplex 0 -dm_refine 0 -temp_petscspace_degree 2 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [2.9]
    suffix: 2d_q2_sconv
    args: -sol_type trig_linear -dm_plex_simplex 0 -temp_petscspace_degree 2 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00000001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 3 -convest_num_refine 3 get L_2 convergence rate: [1.0]
    suffix: 2d_q2_tconv
    args: -sol_type quadratic_trig -dm_plex_simplex 0 -temp_petscspace_degree 2 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu

  test:
    suffix: 3d_p1
    requires: ctetgen
    args: -sol_type quadratic_linear -dm_refine 1 -temp_petscspace_degree 1 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [1.9]
    suffix: 3d_p1_sconv
    requires: ctetgen
    args: -sol_type quadratic_linear -temp_petscspace_degree 1 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 4 -convest_num_refine 3 get L_2 convergence rate: [1.2]
    suffix: 3d_p1_tconv
    requires: ctetgen
    args: -sol_type quadratic_trig -temp_petscspace_degree 1 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    suffix: 3d_p2
    requires: ctetgen
    args: -sol_type quadratic_linear -dm_refine 0 -temp_petscspace_degree 2 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [2.9]
    suffix: 3d_p2_sconv
    requires: ctetgen
    args: -sol_type trig_linear -temp_petscspace_degree 2 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00000001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 3 -convest_num_refine 3 get L_2 convergence rate: [1.0]
    suffix: 3d_p2_tconv
    requires: ctetgen
    args: -sol_type quadratic_trig -temp_petscspace_degree 2 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    suffix: 3d_q1
    args: -sol_type quadratic_linear -dm_plex_simplex 0 -dm_refine 1 -temp_petscspace_degree 1 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [1.9]
    suffix: 3d_q1_sconv
    args: -sol_type quadratic_linear -dm_plex_simplex 0 -temp_petscspace_degree 1 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 4 -convest_num_refine 3 get L_2 convergence rate: [1.2]
    suffix: 3d_q1_tconv
    args: -sol_type quadratic_trig -dm_plex_simplex 0 -temp_petscspace_degree 1 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    suffix: 3d_q2
    args: -sol_type quadratic_linear -dm_plex_simplex 0 -dm_refine 0 -temp_petscspace_degree 2 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 2 -convest_num_refine 3 get L_2 convergence rate: [2.9]
    suffix: 3d_q2_sconv
    args: -sol_type trig_linear -dm_plex_simplex 0 -temp_petscspace_degree 2 -ts_convergence_estimate -ts_convergence_temporal 0 -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 1 -ts_dt 0.00000001 -snes_error_if_not_converged -pc_type lu
  test:
    # -dm_refine 3 -convest_num_refine 3 get L_2 convergence rate: [1.0]
    suffix: 3d_q2_tconv
    args: -sol_type quadratic_trig -dm_plex_simplex 0 -temp_petscspace_degree 2 -ts_convergence_estimate -convest_num_refine 1 \
          -ts_type beuler -ts_max_steps 4 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu

  test:
    # For a nice picture, -bd_dm_refine 2 -dm_refine 1 -dm_view hdf5:${PETSC_DIR}/sol.h5 -ts_monitor_solution hdf5:${PETSC_DIR}/sol.h5::append
    suffix: egads_sphere
    requires: egads ctetgen
    args: -sol_type quadratic_linear \
          -dm_plex_boundary_filename ${wPETSC_DIR}/share/petsc/datafiles/meshes/sphere_example.egadslite -dm_plex_boundary_label marker \
          -temp_petscspace_degree 2 -dmts_check .0001 \
          -ts_type beuler -ts_max_steps 5 -ts_dt 0.1 -snes_error_if_not_converged -pc_type lu

TEST*/
