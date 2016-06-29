AM_CXXFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS) $(LDEXE)

ICONS=icon-cancel.h icon-failed.h icon-monitor.h icon-running.h icon-draw.h	\
icon-finished.h icon-optics.h icon-waiting.h mouse_black.h mouse_hand.h mouse_white.h mygtk.h\
button-err.png button-ok.png button-play.png button-wait.png icon-waiting.png monitor.png mouse_black.png mouse_hand.png mouse_white.png monitor.html
noinst_HEADERS=monitor.h drawdaemon.h scheduler_ws.h $(ICONS) 

PROGS=

if USE_SCHEDULER
PROGS+=../bin/scheduler
___bin_scheduler_SOURCES=scheduler.cpp
___bin_scheduler_LDADD=../sys/libaosys.la $(LDOPTS)
if HAS_LWS
___bin_scheduler_SOURCES+=scheduler_ws.cpp
___bin_scheduler_LDADD+=-lwebsockets
endif
endif

if LINUX
PROGS+=../bin/load ../bin/record_cpu
___bin_load_SOURCES=load.cpp
___bin_record_cpu_SOURCES=record_cpu.cpp
endif

if DRAW
PROGS+=../bin/drawdaemon 
___bin_drawdaemon_SOURCES=drawdaemon.cpp drawdaemon_gui.cpp drawdaemon_io.cpp drawdaemon_draw.cpp
___bin_drawdaemon_CXXFLAGS=$(DRAW_CFLAGS) $(AM_CXXFLAGS)
___bin_drawdaemon_LDADD=../sys/libaosys.la $(DRAW_LIBS) -lgthread-2.0
endif

if USE_SCHEDULER
if MON
if MONN #FULL with notify
MONITOR_FLAGS=$(MONN_CFLAGS) $(AM_CXXFLAGS) -DWITH_NOTIFY=1
MONITOR_LIBS=../sys/libaosys.la $(LDOPTS) $(MONN_LIBS) 
else #Monitor wihtout libnotify.
MONITOR_FLAGS=$(MON_CFLAGS) $(AM_CXXFLAGS) -DWITH_NOTIFY=0
MONITOR_LIBS=../sys/libaosys.la $(LDOPTS) $(MON_LIBS) 
endif
PROGS+=../bin/monitor ../bin/monitor_list ../bin/monitor_table
___bin_monitor_SOURCES=monitor.cpp monitor_gui.cpp monitor_thread.cpp
___bin_monitor_CXXFLAGS=$(MONITOR_FLAGS)
___bin_monitor_LDADD=$(MONITOR_LIBS)
___bin_monitor_list_SOURCES=monitor_gui_list.cpp monitor.cpp monitor_thread.cpp
___bin_monitor_list_CXXFLAGS=$(MONITOR_FLAGS)
___bin_monitor_list_LDADD=$(MONITOR_LIBS)
___bin_monitor_table_SOURCES=monitor_gui_table.cpp monitor.cpp monitor_thread.cpp
___bin_monitor_table_CXXFLAGS=$(MONITOR_FLAGS)
___bin_monitor_table_LDADD=$(MONITOR_LIBS)
endif
endif

LDADD=../lib/libaos.la $(LDEXE_FULL)
PROGS+=../bin/bin2fits ../bin/psfiris ../bin/psfgc ../bin/ethtest ../bin/mvm_cpu ../bin/drawbin ../bin/drawres 
___bin_bin2fits_SOURCES=bin2fits.cpp
___bin_psfiris_SOURCES=psfiris.cpp
___bin_psfgc_SOURCES=psfgc.cpp
___bin_ethtest_SOURCES=ethtest.cpp
___bin_mvm_cpu_SOURCES=mvm_cpu.cpp
___bin_drawbin_SOURCES=drawbin.cpp
___bin_drawres_SOURCES=drawres.cpp
bin_PROGRAMS=$(PROGS)
