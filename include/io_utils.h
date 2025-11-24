// ==================== io_utils.h ====================
#ifndef IO_UTILS_H
#define IO_UTILS_H

#include "Types.h"

void io_write_solution(const DGSolver* solver, const char* filename);

void io_compute_error(const DGSolver* solver);

double io_get_wall_time(void);

void io_write_modal_coefficients(const DGSolver* solver, const char* filename);

void io_write_high_resolution(const DGSolver* solver, const char* filename);

void io_write_diagnostics(const DGSolver* solver, const char* filename, int append);

void io_write_cell_detail(const DGSolver* solver, int cell_id, const char* filename);

void io_print_progress(const DGSolver* solver, double wall_time);



#endif // IO_UTILS_H