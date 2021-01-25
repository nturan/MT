#pragma once

#include "event_generator.h"
#include "parameters.h"
#include <fstream>
#include <iostream>

void EvaluateLoEvents(EventGenerator* eg, Parameters* p);
void EvaluateNloEvents(EventGenerator* eg, Parameters* p);