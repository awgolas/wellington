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
This module contains the atomic configuration and configurations classes.
"""

import abc

from .. import suite as suiteModule
from .. import misc as miscModule
from ..decays import decayData as decayDataModule
from ..quantities import bindingEnergy as bindingEnergyModule

class configuration( miscModule.classWithSubshellKey ) :

    moniker = 'configuration'

    def __init__( self, subshell, electronNumber ) :

        miscModule.classWithSubshellKey.__init__( self, subshell )

        if( not( isinstance( electronNumber, float ) ) ) : raise TypeError( 'electronNumber is not a float' )
        self.__electronNumber = electronNumber

        self.__bindingEnergy = bindingEnergyModule.suite( )
        self.__bindingEnergy.setAncestor( self )

        self.__decayData = decayDataModule.decayData( )
        self.__decayData.setAncestor( self )

    @property
    def electronNumber( self ) :

        return( self.__electronNumber )

    @property
    def bindingEnergy( self ) :

        return( self.__bindingEnergy )

    @property
    def decayData( self ) :

        return( self.__decayData )

    def convertUnits( self, unitMap ) :

        self.__bindingEnergy.convertUnits( unitMap )
        self.__decayData.convertUnits( unitMap )

    def copy( self ) :

        _configuration = self.__class__( self.subshell, self.electronNumber )
        self.__decayData.copyItems( _configuration.decayData )
        return( _configuration )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s subshell="%s" electronNumber="%s">' % ( indent, self.moniker, self.subshell, self.electronNumber ) ]
        XMLStringList += self.__bindingEnergy.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.__decayData.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            if( child.tag == bindingEnergyModule.suite.moniker ) :
                self.bindingEnergy.parseXMLNode( child, xPath, linkData )
            elif( child.tag == decayDataModule.decayData.moniker ) :
                self.decayData.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( 'Invalid tag = "%s"' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( element.attrib['subshell'], float( element.attrib['electronNumber'] ) )
        self.parseXMLNode( element, xPath, linkData )

        xPath.pop( )
        return( self )

class configurations( suiteModule.suite ) :

    moniker = 'configurations'

    def __init__( self, replace = True ) :

        suiteModule.suite.__init__( self, ( configuration, ), replace = replace )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            self.add( configuration.parseXMLNodeAsClass( child, xPath, linkData ) )

        xPath.pop( )
        return( self )
