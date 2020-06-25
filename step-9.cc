/* ---------------------------------------------------------------------
 *
* Copyright (C) 2000 - 2018 by the deal.II authors
*
 * 
  * This file is part of the deal.II example step-9
  *  Solving a advection-diffusion equation  
  * created by SPZ 
 * ---------------------------------------------------------------------
*
 *
  * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/numerics/error_estimator.h>
#include <fstream>
#include <iostream>


namespace Step9
{
  using namespace dealii;
//                                                           Beta 
  template <int dim>
  class AdvectionField : public TensorFunction<1, dim>
  {
  public:
    virtual Tensor<1, dim> value(const Point<dim> &p) const override;

     };

  template <int dim>
  Tensor<1, dim> AdvectionField<dim>::value(const Point<dim> &p) const
  {
    Point<dim> value;
    value[0] = 2;
    for (unsigned int i = 1; i < dim; ++i)
      value[i] = 1 + 0.8 * std::sin(8. * numbers::PI * p[0]);

    return value;
  }
    //                                              F (rhs) 
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;

  private:

  };

 
  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> & p,
                                   const unsigned int component) const
  {
     (void)component;
    Assert(component == 0, ExcIndexRange(component, 0, 1)); 
  double return_value = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value +=  0.01* std::pow(p(i), 2.0);

  return return_value;
  }
 //                                                     Inflow bc
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };

  template <int dim>                                                        
  double BoundaryValues<dim>::value(const Point<dim> & p,
                                    const unsigned int component) const
  {
     (void)component;
    Assert(component == 0, ExcIndexRange(component, 0, 1)); 
   
    const double Vel_inlet    = -std::pow(p(1), 2.0) +1. ;
    return Vel_inlet;
  }

  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem();

    void run();

  private:
    void setup_system();
    
 //  ///                                                    scratch data: temporary storage for paralelization
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<dim> &fe);
      AssemblyScratchData(const AssemblyScratchData &scratch_data);
  
      FEValues<dim>     fe_values;
  
      std::vector<double>         rhs_values;
      std::vector<Tensor<1, dim>> advection_directions;                                       

      // Finally, we need objects that describe the problem's data:
      AdvectionField<dim> advection_field;
      RightHandSide<dim>  right_hand_side;
      BoundaryValues<dim> boundary_values;                          
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    void assemble_system();
    void local_assemble_system(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      AssemblyScratchData &                                 scratch,
      AssemblyCopyData &                                    copy_data);
    void copy_local_to_global(const AssemblyCopyData &copy_data);

    void boundary_condition(); 
    void solve();
	void output_results() const;
	

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
   
    FE_Q<dim> fe;
    
    AffineConstraints<double> hanging_node_constraints;

    SparsityPattern      sparsity_pattern;       // temporary container for sizing and etc.
    SparseMatrix<double> system_matrix;
	
    Vector<double> solution;                    //    A
    Vector<double> system_rhs;                  //    F
  };

 
  template <int dim>                               // Initialization                  
  AdvectionProblem<dim>::AdvectionProblem()
    : dof_handler(triangulation)
    , fe(3)                             // degree of polynomial FEM
  {}



  template <int dim>
  void AdvectionProblem<dim>::setup_system()                   // constituting U and F matrices
  {
    dof_handler.distribute_dofs(fe);
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    hanging_node_constraints,
                                    /*keep_constrained_dofs =*/false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  template <int dim>
  void AdvectionProblem<dim>::assemble_system()
  {
    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &AdvectionProblem::local_assemble_system,
                    &AdvectionProblem::copy_local_to_global,
                    AssemblyScratchData(fe),
                    AssemblyCopyData());
  }


  template <int dim>
  AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<dim> &fe)
    : fe_values(fe,
                QGauss<dim>(fe.degree + 1),
                update_values | update_gradients | update_quadrature_points |
                  update_JxW_values)
    , rhs_values(fe_values.get_quadrature().size())
    , advection_directions(fe_values.get_quadrature().size())                         // beta advec-dire  !

  {}


  template <int dim>
  AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
	
    : fe_values(scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_values | update_gradients | update_quadrature_points |
                  update_JxW_values)
    , rhs_values(scratch_data.rhs_values.size())
    , advection_directions(scratch_data.advection_directions.size())
  {}



  template <int dim>
  void AdvectionProblem<dim>::local_assemble_system(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    AssemblyScratchData &                                 scratch_data,
    AssemblyCopyData &                                    copy_data)
  {
    // We define some abbreviations to avoid unnecessarily long lines:
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points =
      scratch_data.fe_values.get_quadrature().size();
  
	
    // We declare cell matrix and cell right hand side...
    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    copy_data.cell_rhs.reinit(dofs_per_cell);

    // ... an array to hold the global indices of the degrees of freedom of
    // the cell on which we are presently working...
    copy_data.local_dof_indices.resize(dofs_per_cell);

    // ... then initialize the <code>FEValues</code> object...
    scratch_data.fe_values.reinit(cell);

    // ... obtain the values of right hand side and advection directions 
    scratch_data.advection_field.value_list(
      scratch_data.fe_values.get_quadrature_points(),
      scratch_data.advection_directions);
    scratch_data.right_hand_side.value_list(
      scratch_data.fe_values.get_quadrature_points(), scratch_data.rhs_values);

    // ... and assemble the local contributions to the system matrix and
    // right hand side as also discussed above:
    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          // Alias the AssemblyScratchData object to keep the lines from
          // getting too long:
          const auto &sd = scratch_data;
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            copy_data.cell_matrix(i, j) +=                                   // U
              (( sd.fe_values.shape_value(i, q_point) *       // (phi_i *
                 sd.advection_directions[q_point]     *      //  beta .
                 sd.fe_values.shape_grad(j, q_point)) +      // grad phi_j)
                 sd.fe_values.shape_grad(i, q_point)  *      // grad phi_i *
                 sd.fe_values.shape_grad(j, q_point)) *      // grad phi_j
                 sd.fe_values.JxW(q_point);                         // dx

          copy_data.cell_rhs(i) +=                                           // F (right_hand_side)
            sd.fe_values.shape_value(i, q_point) *          // (phi_i *
             sd.rhs_values[q_point] *                       // f *
             sd.fe_values.JxW(q_point);                      // dx
        }

//         indicating boundary id for inflow at left boundary

    for (unsigned int face_n = 0; face_n < GeometryInfo<dim>::faces_per_cell;
         ++face_n)
      if (cell->face(face_n)->at_boundary())  
        {
		
			if (cell->face(face_n)->center()[0] == -1)
	             cell->face(face_n)-> set_boundary_id (1);
        }
		

    cell->get_dof_indices(copy_data.local_dof_indices);


  }
  
 //  ///             inverting local index into global 

  template <int dim>
  void
  AdvectionProblem<dim>::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    hanging_node_constraints.distribute_local_to_global(
      copy_data.cell_matrix,
      copy_data.cell_rhs,
      copy_data.local_dof_indices,
      system_matrix,
      system_rhs);
  }
  

  template <int dim>
  void AdvectionProblem<dim>::solve()
  {
    SolverControl        solver_control(std::max<std::size_t>(1000,
                                             system_rhs.size() / 10),
                                 1e-10 * system_rhs.l2_norm());
    SolverGMRES<>        solver(solver_control);
    PreconditionJacobi<> preconditioner;
    preconditioner.initialize(system_matrix, 1.0);
    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    Vector<double> residual(dof_handler.n_dofs());

    system_matrix.vmult(residual, solution);
    residual -= system_rhs;
    std::cout << "   Iterations required for convergence: "
              << solver_control.last_step() << '\n'
              << "   Max norm of residual:                "
              << residual.linfty_norm() << '\n';

    hanging_node_constraints.distribute(solution);
  }


  template <int dim>
   void AdvectionProblem<dim>::output_results() const       

  {
    {
      GridOut       grid_out;
	  std::ofstream output("grid-2d.vtk");        
      grid_out.write_vtk(triangulation, output);
    }

    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");
      data_out.build_patches(8);

      DataOutBase::VtkFlags vtk_flags;
      vtk_flags.compression_level =
        DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
      data_out.set_flags(vtk_flags);

	  std::ofstream output("solution-2d.vtk");     
      data_out.write_vtk(output);
    }
  }

 //                                 apply inflow boundary condition 
   template <int dim>
   void AdvectionProblem<dim>::boundary_condition() 
   {
	   
    // Left bc boundary_id  :  1 as an inflow 
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           1,
                                           BoundaryValues<dim>(),
                                           boundary_values);						   
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
   }



  template <int dim>
  void AdvectionProblem<dim>::run()
  {
  
	    GridGenerator::hyper_cube(triangulation, -1, 1);
        triangulation.refine_global(3);
	    std::cout << "   Number of active cells:              "
                  << triangulation.n_active_cells() << std::endl;

        setup_system();

        std::cout << "   Number of degrees of freedom:        "
                  << dof_handler.n_dofs() << std::endl;
		
        assemble_system();
        boundary_condition ();
        solve();
	    output_results();
 
  }

} // namespace Step9


int main()
{
  using namespace dealii;

	
      MultithreadInfo::set_thread_limit();
      
      Step9::AdvectionProblem<2> advection_problem_2d;
      advection_problem_2d.run();
	  
	  
  return 0;
}
