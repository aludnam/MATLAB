// keepontop.c, keepontop.mexw32
//
// Usage:    keepontop('full window name')
//
// Desc:     This function will set a window position to "Always On Top" so that
//           clicking on another window will not make it lose visibility. Most
//           useful for wait/progress bars or other small figures.
//
// Author:   Matthew Parent
// Modified: 7/17/07

#include "mex.h"
#include <windows.h>

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
    RECT rect;
    HWND window;
    char *window_name;
    int length;

    length = mxGetN(prhs[0])+1;
    window_name = mxCalloc(length,sizeof(char)); // Preallocate dynamic memory to hold the string

    mxGetString(prhs[0],window_name,length); // Get input string

    window = FindWindow(NULL,window_name); // Find the window with the given name

    if (!IsWindow(window)){ // Error out if FindWindow returns a non-window
        mexErrMsgTxt("Error: Unable to find window with specified name.");}

    ShowWindow(window,SW_SHOWNORMAL);

    GetWindowRect(window, &rect);
    SetWindowPos(window, HWND_TOPMOST, rect.left, rect.top, 0, 0, SWP_NOSIZE);
}
