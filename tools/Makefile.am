AM_CFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS) $(LDEXE)

ICONS=icon-cancel.h icon-failed.h icon-monitor.h icon-running.h icon-draw.h	\
icon-finished.h icon-optics.h icon-waiting.h mouse_black.h mouse_hand.h mouse_white.h mygtk.h\
button-err.png button-ok.png button-play.png button-wait.png icon-waiting.png monitor.png mouse_black.png mouse_hand.png mouse_white.png monitor.html
noinst_HEADERS=monitor.h drawdaemon.h scheduler_ws.h $(ICONS) 

PROGS=

if USE_SCHEDULER
PROGS+=../bin/scheduler
___bin_scheduler_SOURCES=scheduler.c
___bin_scheduler_LDADD=../sys/libaosys.la $(LDOPTS)
if HAS_LWS
___bin_scheduler_SOURCES+=scheduler_ws.c
___bin_scheduler_LDADD+=-lwebsockets
endif
endif

if LINUX
PROGS+=../bin/load ../bin/record_cpu
___bin_load_SOURCES=load.c
___bin_record_cpu_SOURCES=record_cpu.c
endif

if DRAW
PROGS+=../bin/drawdaemon 
___bin_drawdaemon_SOURCES=drawdaemon.c drawdaemon_gui.c drawdaemon_io.c drawdaemon_draw.c
___bin_drawdaemon_CFLAGS=$(DRAW_CFLAGS) $(AM_CFLAGS)
___bin_drawdaemon_LDADD=../sys/libaosys.la $(DRAW_LIBS) -lgthread-2.0
endif

if USE_SCHEDULER
if MON
if MONN #FULL with notify
MONITOR_FLAGS=$(MONN_CFLAGS) $(AM_CFLAGS) -DWITH_NOTIFY=1
MONITOR_LIBS=../sys/libaosys.la $(LDOPTS) $(MONN_LIBS) 
else #Monitor wihtout libnotify.
MONITOR_FLAGS=$(MON_CFLAGS) $(AM_CFLAGS) -DWITH_NOTIFY=0
MONITOR_LIBS=../sys/libaosys.la $(LDOPTS) $(MON_LIBS) 
endif
PROGS+=../bin/monitor 
___bin_monitor_SOURCES=monitor_gui_list.c monitor.c monitor_thread.c
___bin_monitor_CFLAGS=$(MONITOR_FLAGS)
___bin_monitor_LDADD=$(MONITOR_LIBS)
endif
endif

LDADD=../lib/libaos.la $(LDEXE_FULL)
PROGS+=../bin/bin2fits ../bin/psfiris ../bin/psfgc ../bin/ethtest ../bin/mvm_cpu ../bin/drawbin ../bin/drawres 
___bin_bin2fits_SOURCES=bin2fits.c
___bin_psfiris_SOURCES=psfiris.c
___bin_psfgc_SOURCES=psfgc.c
___bin_ethtest_SOURCES=ethtest.c
___bin_mvm_cpu_SOURCES=mvm_cpu.c
___bin_drawbin_SOURCES=drawbin.c
___bin_drawres_SOURCES=drawres.c
bin_PROGRAMS=$(PROGS)
