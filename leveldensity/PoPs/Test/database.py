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

from PoPs import database as databaseModule
from PoPs import alias as aliasModule

from PoPs.groups import misc as chemicalElementMiscModule

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
from PoPs.families import nuclide as nuclideModule

from PoPs.groups import isotope as isotopeModule
from PoPs.groups import chemicalElement as chemicalElementModule

database = databaseModule.database( 'test', '1.2.3' )

def nuclides( Z, A, data ) :

    symbol = chemicalElementMiscModule.symbolFromZ[Z]
    isotopeID = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( symbol, A )
    keys = random.sample( [ key for key in data ], len( data ) )
    for index in keys :
        mass, energy, charge, halflife, spin, parity = data[index]

        name = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeID, index )
        level = nuclideModule.particle( name )
        energy = nuclearEnergyLevelModule.double( 'base', energy, quantityModule.stringToPhysicalUnit( 'keV' ) )
        level.nucleus.energy.add( energy )

        if( mass is not None ) :
            mass = massModule.double( 'base', mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
            level.mass.add( mass )

        if( charge is not None ) :
            charge = chargeModule.integer( 'base', charge, quantityModule.stringToPhysicalUnit( 'e' ) )
            level.charge.add( charge )

        if( halflife is not None ) :
            time, unit = halflife.split( )
            halflife = halflifeModule.double( 'base', float( time ), quantityModule.stringToPhysicalUnit( unit ) )
            level.halflife.add( halflife )

        if( spin is not None ) :
            spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
            level.spin.add( spin )

        if( parity is not None ) :
            parity = parityModule.integer( 'base', parity, quantityModule.stringToPhysicalUnit( '' ) )
            level.parity.add( parity )

        database.add( level )

#                             mass,  energy, charge, halflife, spin, parity
O16Data = { 0 : [ 15.99491461956,       0,   None,     None, None,   None ],
            1 : [           None, 6049.400,   None,     None, None,   None ],
            2 : [           None, 6129.893,   None,     None, None,   None ],
            3 : [           None, 6917.100,   None,     None, None,   None ] }
for A, data in [ [ 17, O16Data ], [ 16, O16Data ] ] : nuclides( 8, A, data )

#                             mass,  energy, charge,  halflife, spin, parity
Am242Data = { 0 : [ 242.059549159,      0,   None, '16.02 h',    1,     -1 ],
              1 : [          None, 44.092,   None,      None,    0,     -1 ],
              2 : [          None,  48.60,   None,  '141 yr',    5,     -1 ],
              3 : [          None,  52.70,   None,      None,    3,     -1 ] }
nuclides( 95, 242, Am242Data )
database.add( aliasModule.metaStable( 'Am242_m1', 'Am242_e2', 1 ) )

photon = gaugeBosonModule.particle( 'photon' )

mass = massModule.double( 'base', 0, quantityModule.stringToPhysicalUnit( 'amu' ) )
photon.mass.add( mass )

charge = chargeModule.integer( 'base', 0, quantityModule.stringToPhysicalUnit( 'e' ) )
photon.charge.add( charge )

halflife = halflifeModule.string( 'base', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
photon.halflife.add( halflife )

spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( '1' ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
photon.spin.add( spin )

parity = parityModule.integer( 'base', 1, quantityModule.stringToPhysicalUnit( '' ) )
photon.parity.add( parity )

database.add( photon )

leptons = [ [ 'e-',      5.48579909070e-4, -1,     'stable', '1/2',  1, 'electronic' ],
            [ 'e-_anti', 5.48579909070e-4,  1,     'stable', '1/2', -1, 'electronic' ],
            [ 'mu',      0.1134289267,     -1, 2.1969811e-6, '1/2',  1, 'muonic' ] ]

for id, mass, charge, halflife, spin, parity, generation in leptons :
    lepton = leptonModule.particle( id, generation = generation )

    mass = massModule.double( 'base', mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
    lepton.mass.add( mass )

    charge = chargeModule.integer( 'base', charge, quantityModule.stringToPhysicalUnit( 'e' ) )
    lepton.charge.add( charge )

    if( halflife == 'stable' ) :
        halflife = halflifeModule.string( 'base', halflife, quantityModule.stringToPhysicalUnit( 's' ) )
    else :
        halflife = halflifeModule.double( 'base', halflife, quantityModule.stringToPhysicalUnit( 's' ) )
    lepton.halflife.add( halflife )

    spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
    lepton.spin.add( spin )

    parity = parityModule.integer( 'base', parity, quantityModule.stringToPhysicalUnit( '' ) )
    lepton.parity.add( parity )

    database.add( lepton )

database.add( aliasModule.particle( 'electron', 'e-' ) )
database.add( aliasModule.particle( 'e+', 'e-_anti' ) )
database.add( aliasModule.particle( 'positron', 'e-_anti' ) )

baryons = [ [ 'n', 1.00866491588,     0,    881.5, '1/2', 1 ],
            [ 'p', 1.007276466812,    1, 'stable', '1/2', 1 ] ]

for _id, _mass, _charge, _halflife, _spin, _parity in baryons :
    for anti in [ '', '_anti' ] :
        baryon = baryonModule.particle( _id + anti )

        mass = massModule.double( 'base', _mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        baryon.mass.add( mass )

        charge = chargeModule.integer( 'base', _charge, quantityModule.stringToPhysicalUnit( 'e' ) )
        baryon.charge.add( charge )

        if( _halflife == 'stable' ) :
            halflife = halflifeModule.string( 'base', _halflife, quantityModule.stringToPhysicalUnit( 's' ) )
        else :
            halflife = halflifeModule.double( 'base', _halflife, quantityModule.stringToPhysicalUnit( 's' ) )
        baryon.halflife.add( halflife )

        spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( _spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
        baryon.spin.add( spin )

        parity = parityModule.integer( 'base', _parity, quantityModule.stringToPhysicalUnit( '' ) )
        baryon.parity.add( parity )

        database.add( baryon )

xmld1 = database.toXML( )
print xmld1
database2 = database.parseXMLStringAsClass( xmld1 )
if( xmld1 != database2.toXML( ) ) : raise Exception( 'Fix me.' )
