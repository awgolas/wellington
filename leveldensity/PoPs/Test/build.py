#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import random
import fractions

from PoPs import database as databaseModule
from PoPs import alias as aliasModule
from PoPs import misc as miscModule
from PoPs import IDs as IDsModule

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import lepton as leptonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nucleus as nucleusModule
from PoPs.families import nuclide as nuclideModule

from PoPs.groups import isotope as isotopeModule
from PoPs.groups import chemicalElement as chemicalElementModule

oneHalf = fractions.Fraction( '1/2' )
one = fractions.Fraction( '1' )
spinUnit = spinModule.baseUnit

database = databaseModule.database( 'test', '1.2.3' )

particle = miscModule.buildParticleFromRawData( leptonModule.particle, IDsModule.electron,
        mass = ( 0.0005485801, 'amu' ),         spin = ( oneHalf, spinUnit ),                     parity = ( 1, '' ),
        charge = ( -1, 'e' ),                   halflife = ( 'stable', 's' ),               generation = 'electronic' )
database.add( particle )

particle = miscModule.buildParticleFromRawData( leptonModule.particle, IDsModule.electron + '_anti',
        mass = ( 0.0005485801, 'amu' ),         spin = ( oneHalf, spinUnit ),                     parity = ( -1, '' ),
        charge = (  1, 'e' ),                   halflife = ( 'stable', 's' ),               generation = 'electronic' )
database.add( particle )

particle = miscModule.buildParticleFromRawData( leptonModule.particle, 'mu',
        mass = ( 0.1134289267, 'amu' ),         spin = ( oneHalf, spinUnit ),                     parity = ( 1, '' ),
        charge = ( -1, 'e' ),                   halflife = ( 2.1969811e-6, 's' ),           generation = 'muonic' )
database.add( particle )

database.add( aliasModule.particle( 'electron', 'e-' ) )
database.add( aliasModule.particle( 'e+', 'e-_anti' ) )
database.add( aliasModule.particle( 'positron', 'e-_anti' ) )

particle = miscModule.buildParticleFromRawData( baryonModule.particle, 'n',
        mass = ( 1.00866491588, 'amu' ),         spin = ( oneHalf, spinUnit ),                    parity = ( 1, '' ),
        charge = (  0, 'e' ),                   halflife = ( 881., 's' ) )
database.add( particle )

particle = miscModule.buildParticleFromRawData( baryonModule.particle, 'p',
        mass = ( 1.007276466812, 'amu' ),         spin = ( oneHalf, spinUnit ),                   parity = ( 1, '' ),
        charge = (  1, 'e' ),                   halflife = ( 'stable', 's' ) )
database.add( particle )

nucleus = miscModule.buildParticleFromRawData( nucleusModule.particle, 'O16', index = 0, energy = ( 0.0, 'eV' ) )
particle = miscModule.buildParticleFromRawData( nuclideModule.particle, 'O16',
        mass = ( 15.994913988, 'amu' ), energy = ( 0.0, 'eV' ) )
database.add( particle )

nucleus = miscModule.buildParticleFromRawData( nucleusModule.particle, 'o16_e3', index = 3, energy = ( 6917100.0, 'eV' ) )
particle = miscModule.buildParticleFromRawData( nuclideModule.particle, 'O16_e3', nucleus = nucleus )
database.add( particle )

xmld1 = database.toXML( )
print xmld1
database2 = database.parseXMLStringAsClass( xmld1 )

if( xmld1 != database2.toXML( ) ) : raise Exception( 'Fix me.' )
