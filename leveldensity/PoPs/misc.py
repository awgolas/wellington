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

from xData import ancestry as ancestryModule

maxLength = 32
antiSuffix = '_anti'

def toLimitedString( object, maxLength = maxLength ) :
    """
    Returns a str instance of object that is truncated to no more than 32 characters.
    """

    string = str( object )
    if( len( string ) <= maxLength ) : return( string )
    return( string[:maxLength-4] + ' ...' )

def baseAntiQualifierFromID( id, qualifierAllowed = False ) :

    if( not( isinstance( id, str ) ) ) : raise TypeError( 'id is not a str: %s' % type( id ) )
    anti, qualifier = '', ''
    base = id.split( '{' )
    if( len( base ) > 1 ) :
        if( len( base ) > 2 ) : raise ValueError( 'Invalid id = "%s"' % toLimitedString( id, 64 ) )
        base, qualifier = base
        if( qualifier[-1] != '}' ) : raise ValueError( 'Invalid qualifier "%s": id = "%s"' % \
                ( toLimitedString( qualifier, 64 ), toLimitedString( id, 64 ) ) )
        qualifier = qualifier[:-1]
    else :
        base = base[0]
    if( base[-5:] == antiSuffix ) : base, anti = base[:-5], base[-5:]

    if( not( qualifierAllowed ) and ( qualifier != '' ) ) : 
        raise ValueError( 'Particle id cannot have a qualifier: id = "%s"' % toLimitedString( id ) )

    return( base, anti, qualifier )

def baseAntiFromID( id ) :

    base, anti, qualifier = baseAntiQualifierFromID( id )

    return( base, anti )

def buildParticleFromRawData( cls, ID, mass = None, spin = None, parity = None, charge = None, halflife = None,
        nucleus = None, index = None, energy = None, generation = None, label = 'default' ) :

    from .quantities import quantity as quantityModule
    from .quantities import charge as chargeModule
    from .quantities import nuclearEnergyLevel as nuclearEnergyLevelModule
    from .families import particle as particleModule
    from .families import lepton as leptonModule
    from .families import nucleus as nucleusModule
    from .families import nuclide as nuclideModule

    def getUnit( unit ) :

        if( isinstance( unit, str ) ) : unit = quantityModule.stringToPhysicalUnit( unit )
        return( unit )

    if( issubclass( cls, leptonModule.particle ) ) :
        particle = cls( ID, generation = generation )
    elif( issubclass( cls, nucleusModule.particle ) ) :
        ID = ID[0].lower( ) + ID[1:]
        if( index is None ) : raise ValueError( 'index must be defined for nuclide to be built' )
        particle = cls( ID, index )
        if( energy is not None ) : particle.energy.add( nuclearEnergyLevelModule.double( label, energy[0], getUnit( energy[1] ) ) )
    elif( issubclass( cls, nuclideModule.particle ) ) :
        particle = cls( ID )
        if( nucleus is not None ) : particle.nucleus.replicate( nucleus )
        if( charge is None ) : charge = ( 0, chargeModule.baseUnit )
        if( len( particle.nucleus.charge ) == 0 ) :
            particle.nucleus.charge.add( chargeModule.integer( label, particle.Z, chargeModule.baseUnit ) )
        if( energy is not None ) : particle.nucleus.energy.add( nuclearEnergyLevelModule.double( label, energy[0], getUnit( energy[1] ) ) )
    elif( issubclass( cls, particleModule.particle ) ) :
        particle = cls( ID )
    else :
        raise TypeError( 'Invalid class.' )

    particle.buildFromRawData( mass = mass, spin = spin, parity = parity, charge = charge, halflife = halflife, label = label )

    return( particle )

def returnAntiParticleIDFromId( particleID ) :

    n1 = len( antiSuffix )
    if( particleID[-n1:] == antiSuffix ) : return( particleID[-n1:] )
    return( particleID + antiSuffix )

def returnAntiParticleID( particle ) :

    return( returnAntiParticleIDFromId( particle.ID ) )

class classWithIDKey( ancestryModule.ancestry ) :

    __keyName = 'id'

    def __init__( self, id ) :

        ancestryModule.ancestry.__init__( self  )

        if( not( isinstance( id, str ) ) ) : raise TypeError( 'id not str' )
        self.__id = id

    @property
    def id( self ) :

        return( self.__id )

    @property
    def key( self ) :

        return( self.__id )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'id must be a string instance.' )
        self.__id = value

    @property
    def keyName( self ) :

        return( self.__keyName )

class classWithSymbolKey( ancestryModule.ancestry ) :

    __keyName = 'symbol'

    def __init__( self, symbol ) :

        ancestryModule.ancestry.__init__( self  )

        if( not( isinstance( symbol, str ) ) ) : raise TypeError( 'symbol not str' )
        self.__symbol = symbol

    @property
    def symbol( self ) :

        return( self.__symbol )

    @property
    def key( self ) :

        return( self.__symbol )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'symbol must be a string instance.' )
        self.__symbol = value

    @property
    def keyName( self ) :

        return( self.__keyName )

class classWithLabelKey( ancestryModule.ancestry ) :

    __keyName = 'label'

    def __init__( self, label ) :

        ancestryModule.ancestry.__init__( self  )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label not str' )
        self.__label = label

    @property
    def label( self ) :
    
        return( self.__label )

    @property
    def key( self ) :
    
        return( self.__label )
        
    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string instance.' )
        self.__label = value

    @property
    def keyName( self ) :

        return( self.__keyName )

class classWithIndexKey( ancestryModule.ancestry ) :

    __keyName = 'index'

    def __init__( self, index ) :

        ancestryModule.ancestry.__init__( self  )

        if( not( isinstance( index, str ) ) ) : raise TypeError( 'index not str' )
        self.key = index

    @property
    def index( self ) :

        return( self.__index )

    @property
    def key( self ) :

        return( self.__index )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'indexmust be a string instance.' )
        self.__index = value

    @property
    def keyName( self ) :

        return( self.__keyName )

class classWithSubshellKey( ancestryModule.ancestry ) :

    __keyName = 'subshell'

    def __init__( self, subshell ) :

        ancestryModule.ancestry.__init__( self  )

        if( not( isinstance( subshell, str ) ) ) : raise TypeError( 'subshell not str' )
        self.__subshell = subshell

    @property
    def subshell( self ) :

        return( self.__subshell )

    @property
    def key( self ) :

        return( self.__subshell )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'subshell must be a string instance.' )
        self.__subshell = value

    @property
    def keyName( self ) :

        return( self.__keyName )
