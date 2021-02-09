#pragma once

#include "event_generator.h"
#include "parameters.h"
#include <fstream>
#include <iostream>

void EvaluateLoEvents(Event e, Parameters* p);
void EvaluateNloEvents(Event e, Parameters* p);