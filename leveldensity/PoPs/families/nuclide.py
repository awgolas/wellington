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

"""
This module contains the nuclear level classes.
"""

from .. import misc as miscModule
from ..quantities import nuclearEnergyLevel as nuclearEnergyLevelModule
from ..groups import misc as chemicalElementMiscModule
from ..fissionFragmentData import fissionFragmentData as fissionFragmentDataModule

from . import particle as particleModule
from . import nucleus as nucleusModule

class alias( particleModule.alias ) :

    moniker = 'nuclideAlias'

    @property
    def chemicalElement( self ) :

        return( self.__particle.chemicalElement )

    @property
    def Z( self ) :

        return( self.__particle.Z )

    @property
    def A( self ) :

        return( self.__particle.A )

    @property
    def index( self ) :

        return( self.__particle.index )

    @property
    def energy( self ) :

        return( self.__particle.energy )

class particle( particleModule.particle ) :

    moniker = 'nuclide'
    alias = alias

    def __init__( self, id ) :

        particleModule.particle.__init__( self, id )

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( id )
        self.__nucleus = nucleusModule.particle( chemicalElementMiscModule.nucleusIDFromNuclideID( id ), levelID )
        self.__nucleus.setAncestor( self, attribute = self.keyName )

        self.__fissionFragmentData = fissionFragmentDataModule.fissionFragmentData( )
        self.__fissionFragmentData.setAncestor( self )

    def __eq__( self, other ) :

        return( self.id == other.id )

    @property
    def fissionFragmentData( self ) :

        return( self.__fissionFragmentData )

    @property
    def nucleus( self ) :

        return( self.__nucleus )

    @property
    def A( self ) :

        return( self.__nucleus.__A )

    @property
    def chemicalElementSymbol( self ) :

        return( self.__nucleus.chemicalElementSymbol )

    @property
    def index( self ) :

        return( self.__nucleus.index )

    @property
    def isotope( self ) :

        return( self.ancestor.ancestor )

    @property
    def energy( self ) :

        return( self.__nucleus.energy )

    @property
    def Z( self ) :

        return( self.__nucleus.Z )

    def check( self, info ):

        from .. import warning as warningModule
        warnings = []

        subWarnings = self.__nucleus.check(info)
        if subWarnings:
            warnings.append( warningModule.context('nucleus', subWarnings) )
        # FIXME other checks to perform on the nuclide? Will decay info ever live in the nuclide?

        return warnings

    def convertUnits( self, unitMap ) :

        particleModule.particle.convertUnits( self, unitMap )
        self.__nucleus.convertUnits( unitMap )
        self.__fissionFragmentData.convertUnits( unitMap )

    def copy( self ) :

        _particle = particle( self.id )
        self.__copyStandardQuantities( _particle )
        _particle.__nucleus.replicate( self.__nucleus )
        _particle.__fissionFragmentData.replicate( self.__fissionFragmentData )

        return( _particle )

    def extraXMLElements( self, indent, **kwargs ) :

        XMLStringList = self.__nucleus.toXMLList( indent, **kwargs )
        XMLStringList += self.__fissionFragmentData.toXMLList( indent, **kwargs )

        return( XMLStringList )

    def getMass( self, unit ) :

        if( len( self.mass ) > 0 ) : return( self.mass[0].float( unit ) )
        if( self.index == 0 ) : raise Exception( 'Recursion detected as group-state does not have a mass: ID = %s.' % self.id )
        return( self.ancestor[0].mass[0].float( unit ) + self.energy[0].float( unit + ' * c**2' ) )

    def parseExtraXMLElement( self, element, xPath, linkData ) :

        if( element.tag == nucleusModule.particle.moniker ) :
            nucleus = nucleusModule.particle.parseXMLNodeAsClass( element, xPath, linkData )
            self.__nucleus.replicate( nucleus )
            return( True )
        elif( element.tag == fissionFragmentDataModule.fissionFragmentData.moniker ) :
            fissionFragmentData = fissionFragmentDataModule.fissionFragmentData.parseXMLNode( element, xPath, linkData )
            self.__fissionFragmentData.replicate( fissionFragmentData )
            return( True )

        return( False )

    def sortCompare( self, other ) :

        if( not( isinstance( other, particle ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.index - other.index )

class suite( particleModule.suite ) :

    moniker = 'nuclides'
    particle = particle
