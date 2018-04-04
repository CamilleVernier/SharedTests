#ifndef SETTINGS_H
#define SETTINGS_H

#include "lattice_s.h"

int read_settings_file(const char filename[]);
int apply_final_settings();
int cmp_nocase_no__(const std::string& s, const std::string& s2);

#endif
