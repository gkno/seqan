//
// $Id: clock.inl,v 1.2 2006/06/17 13:18:45 petzold Exp $
//
// $Log: clock.inl,v $
// Revision 1.2  2006/06/17 13:18:45  petzold
// added const to some function members of class Clock. Is there a better way to get along with constant references\?\!
//
// Revision 1.1  2006/03/09 10:28:59  pmay
// new theseus files
//
// Revision 1.5  2004/04/16 19:10:38  bzcstein
// precision in clock output stream
//
// Revision 1.4  2004/04/09 15:23:43  bzcstein
// get clock ticks from sysconf
//
// Revision 1.3  2003/11/07 09:57:06  bzcstein
// def is wall, usr & sys output
//
// Revision 1.2  2003/10/10 21:10:54  bzcstein
// correct if for wall, usr/sys methods added
//
//

#ifndef _WIN32
#include <unistd.h>
#include <sys/times.h>
#endif

inline Clock::Clock() {
	#ifndef _WIN32
		_clk_tcks = sysconf(_SC_CLK_TCK);
	#endif
	Reset();
}

inline void Clock::Reset() {
#ifdef _WIN32
    this->_tot_time=0;
#else 
    _tot_time = 0;
    _tot_utime = 0;
    _tot_stime = 0;
#endif
}

inline void Clock::Start() {
#ifdef _WIN32
    this->_start=clock();
#else 
    tms my_tms;
    _start  = times(&my_tms);
    _ustart = my_tms.tms_utime + my_tms.tms_cutime;
    _sstart = my_tms.tms_stime + my_tms.tms_cstime;
#endif
}

inline void Clock::Stop() {
#ifdef _WIN32
    this->_end=clock();
	this->_tot_time += this->_end - this->_start;
#else 
    tms my_tms;
    _end  = times(&my_tms);
	_tot_time += _end - _start;
    _uend = my_tms.tms_utime + my_tms.tms_cutime;
	_tot_utime += _uend - _ustart;
    _send = my_tms.tms_stime + my_tms.tms_cstime;
	_tot_stime += _send - _start;
#endif
}


inline std::ostream& operator<<(std::ostream& os,Clock& c) {
    os << std::setiosflags ( std::ios::fixed );
    os << std::setprecision(2);
#ifdef _WIN32
    os << (double)(c._tot_time)/(double)CLOCKS_PER_SEC << "s";
#else 
    os  //<< " wall: " << std::setw(10) << (double)(c._tot_time)  /(double)c._clk_tcks << " s"
        // << " usr: " 
		<< std::setw(10) << (double)( c._tot_utime )/(double)c._clk_tcks << " s";
        //<< " sys: " << std::setw(10) << (double)(c._tot_stime)/(double)c._clk_tcks << " s";
#endif		/* ! _WIN32 */
    os << std::resetiosflags ( std::ios::fixed );
    os << std::setprecision(6);
    return os;
}


inline double Clock::Wall() const {
#ifdef _WIN32
	return ( (double)(this->_tot_time)/(double)CLOCKS_PER_SEC );
#else
	return ( (double)(this->_tot_stime)/(double)_clk_tcks );
#endif
}

inline double Clock::Usr() const {
#ifdef _WIN32
	return ( (double)(this->_tot_time)/(double)CLOCKS_PER_SEC );
#else
	return ( (double)(this->_tot_utime)/(double)_clk_tcks );
#endif
}

inline double Clock::Sys() const {
#ifdef _WIN32
	return ( (double)(this->_tot_time)/(double)CLOCKS_PER_SEC );
#else
	return ( (double)(this->_tot_stime)/(double)_clk_tcks );
#endif
}


