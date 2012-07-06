/*
Copyright (c) 2009, Timothy A. Davis
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

* Redistributions of source code must retain the above copyright 
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in 
the documentation and/or other materials provided with the distribution
* Neither the name of the University of Florida nor the names 
of its contributors may be used to endorse or promote products derived 
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
					   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
					   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
 */

#include "WaitMex.h"

/* -------------------------------------------------------------------------- */
/* waitbar_create: create a waitbar */
/* -------------------------------------------------------------------------- */

/* Create a waitbar and return a point to the new waitbar.  Example:
 * h = waitbar_create (x, "Please wait...") in C is just like
 * h = waitbar (x, 'Please wait') in MATLAB, where x is a fraction between
 * 0 and 1.
 */

waitbar *waitbar_create         /* return a pointer to the new waitbar */
(
 double fraction,            /* fraction from 0 to 1 */
 char *message               /* message to display */
 )
{
    int error ;
    waitbar *h ;
	
    h = mxMalloc (sizeof (waitbar)) ;
    h->fraction = mxCreateDoubleScalar (fraction) ;
    h->message = mxCreateString (message) ;
	
    /* h = waitbar (fraction, message) ; */
    h->inputs [0] = h->fraction ;
    h->inputs [1] = h->message ;
    error = mexCallMATLAB (1, h->outputs, 2, h->inputs, "waitbar") ;
    if (error)
    {
        mexErrMsgTxt ("unable to create waitbar") ;
    }
	
    /* save the MATLAB handle h in the waitbar struct */
    h->handle = h->outputs [0] ;
	
    return (h) ;
}

/* -------------------------------------------------------------------------- */
/* waitbar_update: update a waitbar */
/* -------------------------------------------------------------------------- */

/* Update the length of the bar in an existing waitbar.  Example:
 * waitbar_update (x, h, NULL) in C is just like waitbar (x,h) in MATLAB,
 * where h is the handle to the existing waitbar.  The message is not changed.
 * To change the message, use waitbar_update (x, h, "new message"), which is
 * just like waitbar (x, h, 'new message') in MATLAB.
 */

void waitbar_update
(
 double fraction,
 waitbar *h,
 char *message
 )
{
    int error ;
    if (h == NULL) return ;                 /* nothing to do */    
    (* mxGetPr (h->fraction)) = fraction ;  /* update the fraction */
    h->inputs [0] = h->fraction ;           /* define the inputs x and h */
    h->inputs [1] = h->handle ;
	
    if (message == NULL)
    {
        /* use the existing message; waitbar (x,h) in MATLAB */
        error = mexCallMATLAB (0, h->outputs, 2, h->inputs, "waitbar") ;
    }
    else
    {
        /* define a new message; waitbar (x,h,message) in MATLAB */
        mxDestroyArray (h->message) ;
        h->message = mxCreateString (message) ;
        h->inputs [2] = h->message ;
        error = mexCallMATLAB (0, h->outputs, 3, h->inputs, "waitbar") ;
    }
    if (error)
    {
        mexErrMsgTxt ("unable to update waitbar") ;
    }
}

/* -------------------------------------------------------------------------- */
/* waitbar_destroy: destroy a waitbar */
/* -------------------------------------------------------------------------- */

/* Destroys a waitbar; same as close(h) in MATLAB */

void waitbar_destroy
(
 waitbar *h
 )
{
    int error ;
    if (h == NULL) return ;             /* nothing to do */    
    h->inputs [0] = h->handle ;
    mxDestroyArray (h->fraction) ;      /* free the internal mxArrays */
    mxDestroyArray (h->message) ;
    error = mexCallMATLAB (0, h->outputs, 1, h->inputs, "close") ;
    mxDestroyArray (h->handle) ;
    mxFree (h) ;
    if (error)
    {
        mexErrMsgTxt ("error closing waitbar") ;
    }
}

/* -------------------------------------------------------------------------- */
/* waitbar_return: return a waitbar handle to the caller */
/* -------------------------------------------------------------------------- */

/* This function frees the space used internally in a mexFunction for managing
 * the waitbar, and returns the mxArray handle to the caller.  The waitbar still
 * exists in MATLAB.  Example: pargout [0] = waitbar_return (h) to return the
 * MATLAB handle to the caller of your mexFunction.
 */

mxArray *waitbar_return
(
 waitbar *h
 )
{
    mxArray *handle ;
    if (h == NULL) return (NULL) ;      /* nothing to do */    
    handle = h->handle ;                /* get the MATLAB handle */
    mxDestroyArray (h->fraction) ;      /* free the internal mxArrays */
    mxDestroyArray (h->message) ;
    mxFree (h) ;
    return (handle) ;                   /* return the MATLAB handle */
}
