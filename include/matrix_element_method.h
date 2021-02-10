#pragma once

#include "event_generator.h"
#include "parameters.h"
#include <fstream>
#include <iostream>

void EvaluateLoEvents(std::vector<Event> events, Parameters* p);
void EvaluateNloEvents(std::vector<Event> events, Parameters* p);