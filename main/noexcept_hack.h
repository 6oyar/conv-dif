#ifndef no_except_work_around_H
#	define no_except_work_around_H

#	ifdef _MSC_VER 
#		if (_MSC_VER <= 1800)
#			include <xkeycheck.h>
#			define noexcept throw()
#		endif
#	endif

#endif //no_except_work_around_H
