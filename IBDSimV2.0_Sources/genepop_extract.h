/***************************************************************************
© R. Leblois 2001-2004
© F. Rousset 2005-2008

rousset@isem.univ-montp2.fr

This file is part of IBDSim. This software is a computer program
whose purpose is to perform population genetic simulations.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/
#ifndef H_GENEPOP
#define H_GENEPOP
#include <string>
#include <cstdlib>
#include <cstdio>

extern bool inputCheckBool,batchDebug;
extern std::string EOLtype;
extern const std::string fichierIn;

void _gotoxy(int x,int y);
void effacer_ecran();
int cmp_nocase(const std::string& s, const std::string& s2);

// CBR: previously known as CFichier_genepop::set_eof_check_EOLtype(bool set_eof /*=true*/)
int set_eof_check_EOLtype(std::string EOLFileName, bool set_eof);
// CBR: previously part of CFichier_genepop::parseFile()
void set_UNIX_EOLtype(std::string EOLFileName);

#endif
