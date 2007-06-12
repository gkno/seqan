//
// $Id: clock.h,v 1.3 2006/06/17 13:18:45 petzold Exp $
//
// $Log: clock.h,v $
// Revision 1.3  2006/06/17 13:18:45  petzold
// added const to some function members of class Clock. Is there a better way to get along with constant references\?\!
//
// Revision 1.2  2006/04/03 12:38:11  gunnar
// * changed handling of strstream error (strstream.h now in misc)
// * adapted Makefiles
// * some RAPTOR changes (minor)
// * other minor stuff
//
// Revision 1.1  2006/03/09 10:28:59  pmay
// new theseus files
//
// Revision 1.4  2005/02/07 14:39:26  bzcmaypa
// 64 bit version: Long
//
// Revision 1.3  2004/04/09 15:24:07  bzcstein
// clock ticks from sysconf
//
// Revision 1.2  2003/10/10 21:10:13  bzcstein
// correct wall timing, addition usr/sys methods
//
//

#ifndef CLOCK_INCLUDED
#define CLOCK_INCLUDED
//#include <unistd.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include "config.h"

class Clock {
protected:
    clock_t                 _start;
    clock_t                 _end;
	clock_t                 _tot_time;
    clock_t                 _ustart;
    clock_t                 _uend;
	clock_t                 _tot_utime;
    clock_t                 _sstart;
    clock_t                 _send;
	clock_t                 _tot_stime;
	Long                    _clk_tcks;
    
public:
                            Clock();
	void                    Stop();
    void                    Start();
    void                    Reset();
    friend std::ostream&    operator<<(std::ostream&,Clock&);
    double                  Wall() const;
    double                  Usr() const;
    double                  Sys() const;
};

#include "clock.inl"
#endif

