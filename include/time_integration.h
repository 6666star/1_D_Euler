// ==================== time_integration.h ====================
#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include "Types.h"

void rk3_init_time(TimeInfo* time, double t_end, double dt);
void rk3_advance_time(DGSolver* solver);

#endif // TIME_INTEGRATION_H