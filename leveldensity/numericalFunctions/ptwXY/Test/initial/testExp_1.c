/*
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
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>

#include <nfut_utilities.h>
#include <ptwXY.h>

#define size 1001

static int verbose = 0;

double getDouble( char const *s );
ptwXYPoints *randomUV( statusMessageReporting *smr, double accuracy, int biSectionMax );
void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int doRandom = 0, iarg, echo = 0, biSectionMax = 10;
    int64_t i, j, n;
    ptwXYPoints *u, *y, *e, *p;
    double xyPoints[2*2], x, dx, z, accuracy = 1e-3;
    ptwXYPoint *xy1, *xy2;
    FILE *ff;
    char fmt[] = "%.14e %.14e\n";
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-r", argv[iarg] ) == 0 ) {
            doRandom = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );
    
    nfu_setMemoryDebugMode( 0 );

    xyPoints[0] = 0.0;
    xyPoints[1] = 1.0;
    xyPoints[2] = 1.0;
    xyPoints[3] = -0.2;
    if( argc == 5 ) {
        xyPoints[0] = getDouble( argv[1] );
        xyPoints[1] = getDouble( argv[2] );
        xyPoints[2] = getDouble( argv[3] );
        xyPoints[3] = getDouble( argv[4] );
    }

    ff = fopen( "info.dat", "w" );
    fprintf( ff, "# accuracy = %e\n", accuracy );
    fprintf( ff, "# biSectionMax = %d\n", biSectionMax );
    fclose( ff );

    if( doRandom ) {
        u = randomUV( &smr, accuracy, biSectionMax ); }
    else {
        if( ( u = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, 2, xyPoints, 0 ) ) == NULL ) 
            nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    if( ( y = ptwXY_clone( &smr, u ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );

    ff = fopen( "curve_u.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) u->length );
    ptwXY_simpleWrite( u, ff, fmt );
    fclose( ff );

    if( ptwXY_exp( &smr, y, 1 ) != nfu_Okay ) nfut_printSMRError2p( &smr, "Via." );
    if( ptwXY_simpleCoalescePoints( &smr, y ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    ff = fopen( "exp_u.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) y->length );
    ptwXY_simpleWrite( y, ff, fmt );
    fclose( ff );

    if( ( e = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    n = ptwXY_length( &smr, u );
    xy1 = ptwXY_getPointAtIndex_Unsafely( u, 0 );
    for( i = 1; i < n; i++ ) {
        xy2 = ptwXY_getPointAtIndex_Unsafely( u, i );
        x = xy1->x;
        dx = ( xy2->x - xy1->x ) / 20;
        for( j = 0; j < 20; j++ ) {
            if( ptwXY_getValueAtX( &smr, u, x, &z ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
            if( ptwXY_setValueAtX( &smr, e, x, exp( z ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
            x += dx;
        }
        xy1 = xy2;
    }
    if( ptwXY_getValueAtX( &smr, u, xy2->x, &z ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_setValueAtX( &smr, e, x, exp( z ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    ff = fopen( "exactExp_u.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) e->length );
    ptwXY_simpleWrite( e, ff, fmt );
    fclose( ff );

    if( ptwXY_exp( &smr, u, -1 ) != nfu_Okay ) nfut_printSMRError2p( &smr, "Via." );
    if( ( p = ptwXY_mul2_ptwXY( &smr, u, y ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via. exp( y ) * exp( -y )" );

    ff = fopen( "exp_u_times_exp_minusU.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) p->length );
    ptwXY_simpleWrite( p, ff, fmt );
    fclose( ff );

    ptwXY_free( y );
    ptwXY_free( u );
    ptwXY_free( e );
    ptwXY_free( p );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
double getDouble( char const *s ) {

    double d;
    char *e;

    errno = 0;
    d = strtod( s, &e );
    if( ( *e != 0 ) || ( errno != 0 ) ) printMsg( "could not convert '%s' to double, err = %d, e = %s", s, errno, e );
    return( d );
}
/*
****************************************************************
*/
ptwXYPoints *randomUV( statusMessageReporting *smr, double accuracy, int biSectionMax ) {

    int64_t i;
    double x, y;
    ptwXYPoints *f;

    if( ( f = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    for( i = 0, x = 0, y = 0; i < size; i++ ) {
        x += drand48( );
        y += drand48( ) - 0.5;
        if( ptwXY_setValueAtX( smr, f, x, y ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    }
    return( f );
}
/*
****************************************************************
*/
void printMsg( char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
