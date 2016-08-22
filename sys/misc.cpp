/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <utime.h>
#include <fcntl.h>           /* For O_* constants */
#include <limits.h>
#include <errno.h>
#include <dirent.h>
#include <ctype.h>
#include <stdint.h>
#include "common.h"
#include "thread.h"
#include "process.h"
#include "misc.h"
#include "path.h"
#include "bin.h"
/**
   Obtain the basename of a file. The returnned string must be freed.
*/
char *mybasename(const char *fn){
    if(!fn || strlen(fn)==0) return NULL;
    char fn2[PATH_MAX];
    strncpy(fn2,fn, PATH_MAX);
    /*If this is a folder, remove the last / */
    if(fn2[strlen(fn2)-1]=='/')
	fn2[strlen(fn2)-1]='\0';
    char* sep=strrchr(fn2,'/');
    if(!sep){
	sep=fn2; 
    }else{
	sep++;
    }
    char *bn=(char*)malloc(strlen(sep)+1);
    strcpy(bn,sep);
    return bn;
}
/**
   Copy a file from file stream src to dest.
*/
static void copyfile_fp(FILE *dest, FILE *src){
    char buffer[4096];
    size_t br,bw;
    while(!feof(src)){
	br=fread(buffer, 1, 4096, src);
	if((bw=fwrite(buffer, 1, br, dest))!=br){
	    error("copyfile: Write failed %ld of %ld written.\n", 
		  (long)bw,(long)br);
	}
    }
}
/**
   Copy a file from src to dest
*/
void copyfile(const char *dest, const char *src){
    FILE *psrc=fopen(src,"rb");
    if(!psrc){
	error("Open source file failed\n");
    }
    FILE *pdest=fopen(dest,"wb");
    if(!pdest){
	error("Open destination file failed\n");
    }
    copyfile_fp(pdest, psrc);
}


/**
   Check the suffix of a file.
*/
int check_suffix(const char *fn, const char *suffix){
    if(!fn || !suffix) return 0;
    int lfn=strlen(fn);
    int lsu=strlen(suffix);
    if(lfn < lsu) return 0;
    if(mystrcmp(fn+lfn-lsu,suffix)){
	return 0;
    }else{
	return 1;
    }
}
/**
   Convert argc, argv to a single string, prefixed by the current directory.
*/
char *argv2str(int argc, const char *argv[], const char* delim){
    if(!argc) return NULL;
    char *cwd=mygetcwd();
    size_t slen=strlen(cwd)+2+strlen(HOME);
    if(!delim) delim=" ";
    for(int iarg=0; iarg<argc; iarg++){
	slen+=strlen(delim)+strlen(argv[iarg]);
    }
    char *scmd=mycalloc(slen,char);
    if(!mystrcmp(cwd,HOME)){
	strcpy(scmd,"~");
	strcat(scmd,cwd+strlen(HOME));
    }else{
	strncpy(scmd,cwd, slen);
    }
    strcat(scmd,"/");
    char *exename=mybasename(argv[0]);
    strcat(scmd, exename);
    strcat(scmd, delim);
    free(exename);
    for(int iarg=1; iarg<argc; iarg++){
	if(argv[iarg] && strlen(argv[iarg])>0){
	    strcat(scmd, argv[iarg]);
	    strcat(scmd, delim);
	}
    }
    if(strlen(scmd)>slen-1) error("Overflow\n");
    free(cwd);
    char *scmdend=scmd+strlen(scmd);
    if(delim[0]!='\n'){
	for(char *p=scmd; p<scmdend; p++){
	    if(p[0]=='\n') p[0]=' ';
	}
    }
    while(scmdend>scmd+1 && scmdend[0]==delim[0] && scmdend[-1]==delim[0]){
	scmdend[0]='\0';
	scmdend--;
    }
    return scmd;
}
/**
   Print the content of a file.
*/
void print_file(const char *fnin){
    char *fn=search_file(fnin);
    if(!fn){
	warning("%s not found\n", fnin);
	return;
    }
    FILE *fp;
    if(!(fp=fopen(fn,"r"))){
	error("Open %s failed\n",fn);
    }
    copyfile_fp(stderr, fp);
    fflush(stderr);
    fclose(fp);
    free(fn);
}

/**
   Get current time in seconds as an integer.
*/
int myclocki(){
    time_t a;
    return time(&a);
}
/**
   Get current time in ascii string for easy print. The string
   contains spaces and is not suitable to use in filename. The
   returned string should not be modified.  */
const char *myasctime(void){
    static char st[64];
    time_t a;
    time(&a);
    ctime_r(&a, st);
    st[strlen(st)-1]='\0';/*remove final \n */
    return st;
}
/**
   Get furrent time in ascii string that doesn't contain
   spaces. Suitable for use in filenames. The returnned string
   must be freed. */
char *strtime(void){
    char str[64];
    time_t t=myclocki();
    struct tm tmp;
    localtime_r(&t,&tmp);/*don't free tmp */
    strftime(str,64,"%F-%H%M%S",&tmp);
    char pid[20];
    snprintf(pid, 20, "-%05d", (int)getpid());
    char *dir=stradd(str, pid, NULL);
    return dir;
}

/**
   Get current time in milli-second resolution.
*/
double myclockd(void){
#if defined(_POSIX_TIMERS) && _POSIX_TIMERS > 0
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec+(double)t.tv_nsec*1e-9;
#else
    struct timeval tk;
    gettimeofday(&tk,NULL);
    return (double)tk.tv_sec+(double)tk.tv_usec*1e-6;
#endif
}
/**
   Get current directory. The returnned string must be freed.
*/
char *mygetcwd(void){
    char cwd0[PATH_MAX];
    if(!getcwd(cwd0,PATH_MAX)) 
	error("Error getting current directory\n");
    return strdup(cwd0);
}
/**
   Translate a path into absolute path. The caller shall free the returned string.
*/
char *myabspath(const char *path){
    if(!path) return 0;
    char path2[PATH_MAX];
    switch(path[0]){
    case '~':
	snprintf(path2, PATH_MAX, "%s/%s", HOME, path+1);
	break;
    case '/':
	snprintf(path2, PATH_MAX, "%s", path);
	break;
    default:
	{
	    char *cwd=mygetcwd();
	    snprintf(path2, PATH_MAX, "%s/%s", cwd, path);
	    free(cwd);
	}
    }
    if(!isdir(path)){
	error("path %s does not exist\n", path);
	return 0;
    }else{
	return strdup(path2);
    }
}

void mysymlink(const char *fn, const char *fnlink){
    if(!exist(fn)) return;
    remove(fnlink);
    if(symlink(fn, fnlink)){
	warning("Unable to make symlink %s->%s\n",fnlink,fn);
    }
}
/**
   Test whether a file exists.
*/
int exist(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf);
}

/**
   Test whether fn is directory
*/
int isdir(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISDIR(buf.st_mode);
}

/**
   Test whether fn is ordinary file
*/
int isfile(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISREG(buf.st_mode);
}

/**
   Test whether fn is a symbolic link
*/
int islink(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return !stat(fn, &buf) && S_ISLNK(buf.st_mode);
}
/**
   Test whether fd is a socket
*/
int issock(int fd){
    if(fd==-1) return 0;
    struct stat buf;
    return !fstat(fd, &buf) && S_ISSOCK(buf.st_mode);
}
/**
 * Compute length of file in Bytes
 */
size_t flen(const char *fn){
    if(!fn) return 0;
    struct stat buf;
    return stat(fn, &buf)?0:buf.st_size;
}
/**
   Return the modification time of the file
 */
time_t fmtime(const char *fn){
    struct stat buf;
    if(!fn || stat(fn, &buf)) return 0;
    return buf.st_ctime;
}
/**
   Concatenate many strings. Argument list must end with NULL.
*/
char *stradd(const char* a, ...){
    char *out;
    int n=0;
    va_list ap;
    va_start(ap,a);
    if(a){
	n=strlen(a)+1;
    }
    for(const char *arg=va_arg(ap,const char*); arg; arg=va_arg(ap,const char*)){
	n+=strlen(arg);
    }
    va_end(ap);
    out=mycalloc(n,char);
    if(a){
	strcpy(out,a);
    }
    va_start(ap,a);
    for(const char *arg=va_arg(ap,const char*); arg; arg=va_arg(ap,const char*)){
	strcat(out,arg);
    }
    va_end(ap);
    return out;
}
/**
   Concatenate many strings, like stradd, but arguments are an array of char*
*/
char *strnadd(int argc, const char *argv[], const char* delim){
    if(!argc) return NULL;
    int slen=1+strlen(delim);
    for(int iarg=0; iarg<argc; iarg++){
	slen+=strlen(delim)+strlen(argv[iarg]);
    }
    char *scmd=mycalloc(slen,char);
    for(int iarg=0; iarg<argc; iarg++){
	if(argv[iarg] && strlen(argv[iarg])>0){
	    strcat(scmd,argv[iarg]);
	    strcat(scmd,delim);
	}
    }
    char *scmdend=scmd+strlen(scmd)-1;
    while(scmdend>scmd+1 && scmdend[0]==delim[0] && scmdend[-1]==delim[0]){
	scmdend[0]='\0';
	scmdend--;
    }
    return scmd;
}
/**
   translate a filename into absolute file name that starts with /
*/
char *expand_filename(const char *fn){
    char *out;
    if(fn[0]=='~'){
	out=stradd(HOME,fn+1,NULL);
    }else{
	out=strdup(fn);
    }
    return out;
}
/**
   Duplicate a string. Check for NULL.
 */
char *mystrndup(const char *A, int len){
    int len2=strlen(A);
    if(len2<len) len=len2;
    char *B=(char*)malloc(len+1);
    memcpy(B,A,len);
    B[len]='\0';
    return B;
}

/**
   declare strdup so my memory mangement mem.c is happy when DEBUG=1. Handles
NULL pointer correctly.  */
char *mystrdup(const char *A){
    if(!A){
	return NULL;
    }else{
	int nlen=strlen(A);
	char *B=(char*)malloc(nlen+1);
	memcpy(B,A,nlen+1);
	return B;
    }
}


/**
   Remove files that are older than sec seconds in folder fndir. If sec==0,
   remove everything.
*/
void remove_file_older(const char *fndir, long sec){
    DIR *dir=opendir(fndir);
    if(!dir){
	error("Unable to open directory %s\n",fndir);
    }
    struct dirent *dp;
    struct stat buf;
    char fnfull[PATH_MAX];
    long sec2=myclocki()-sec;
    while((dp=readdir(dir))){
	snprintf(fnfull,PATH_MAX,"%s/%s",fndir,dp->d_name);
	if(!stat(fnfull,&buf) && S_ISREG(buf.st_mode) && (buf.st_atime<=sec2 || sec==0)){
	    remove(fnfull);
	    info2("Remove %s. %ld days old\n", fnfull, (long)(myclocki()-buf.st_mtime)/3600/24);
	}
    }
    closedir(dir);
}

/**
   Make dirs recursively. like mkdir -p in bash
*/
void mymkdir(const char *format, ...){
    format2fn;
    if(!fn) return;
    if(fn[strlen(fn)-1]=='/')
	fn[strlen(fn)-1]='\0';
    if(mkdir(fn, 0700)==-1 && errno!=EEXIST){
	perror("mkdir");
	char *tmp=strrchr(fn,'/');
	if(!tmp){
	    error("Unable to mkdir '%s'\n",fn);
	}
	tmp[0]='\0';
	mymkdir("%s",fn);
	tmp[0]='/';
	if(mkdir(fn,0700)==-1 && errno!=EEXIST){
	    error("Unable to mkdir '%s'\n",fn);
	}
    }
}
/**
   Compare two strings upto the length of b. if length of a is less than b,
   return false. 1 means not equal.
 */
int mystrcmp(const char *a, const char *b){
    if(!a || !b) return 1;
    int la=strlen(a);
    int lb=strlen(b);
    if(la==0 || lb==0 || la<lb){
	return 1;
    }else{
	return strncmp(a,b,lb);
    }
}

/**
   Make the fd close on exec.
*/
void cloexec(int fd){
    int oldflag=fcntl(fd, F_GETFD, 0);
    if(oldflag!=-1){
	oldflag |= FD_CLOEXEC;
	fcntl(fd, F_SETFD, oldflag);
    }
}
/**
   wrap of nanosleep
*/
void mysleep(double sec){
    struct timespec ts;
    ts.tv_sec=(time_t)trunc(sec);
    ts.tv_nsec=(long)((sec-ts.tv_sec)*1e9);
    nanosleep(&ts, NULL);
}

/**
   Return available space of mounted file system in bytes.
*/
long available_space(const char *path){
    struct statvfs buf;
    if(statvfs(path, &buf)){
	perror("statvfs");
	return 0;
    }else{
	return (long)buf.f_bsize * (long)buf.f_bavail;
    }
}
/**
   Extract a string constant from the command line, and output the position
   where the string terminates.*/
static char *cmd_string(char *start, char **end2){
    char *input=start;
    char *end;
    while(isspace((int)input[0]) || input[0]=='\n') input++;
    if(input[0]=='\'' || input[0]== '"'){
	end=strchr(input+1, input[0]);/*find matching quote. */
	input[0]=' ';
	input++;
	if(!end){
	    error("String does not end\n");
	}
    }else{
	end=input;
	while(!isspace((int)end[0]) && end[0]!='\n' && end[0]!='\0') end++;
    }
    int noteos=0;
    if(end[0]!='\0'){
	end[0]='\0';
	noteos=1;
    }
    char *out=strdup(input);
    memset(input, ' ', strlen(input));
    if(noteos){
	end[0]=' ';
	*end2=end+1;
    }else{
	*end2=end;
    }
    return out;
}
/**
   Parse command line arguments. The remaining string contains whatever is not yet parsed. 
   This is more relaxed than the built in getopd
*/
void parse_argopt(char *cmds, ARGOPT_T *options){
    char *cmds_end=cmds+(cmds?strlen(cmds):0);
    char *start=cmds;
    while(start<cmds_end){
        if(isspace((int)start[0]) || start[0]=='\n'){
	    start[0]=' ';
	    start++;
	    continue;
	}
	if(options && start[0]=='-'){
	    char *start0=start;
	    char key='0';
	    char *value;
	    int iopt=-1;
	    start++;
	    if(start[0]=='-'){/*long option, replace with short ones. */
		start++;
		for(int i=0; (options[i].name); i++){
		    if(!mystrcmp(start, options[i].name)){
			key=options[i].key;
			start+=strlen(options[i].name);
			while(isspace((int)start[0]) || start[0]=='\n'){
			    start[0]=' ';
			    start++;
			}
			if(start[0]=='=') {
			    start[0]=' ';
			    start++;
			}
			iopt=i;
			break;
		    }
		}
	    }else{
		key=start[0];
		start++;
		for(int i=0; (options[i].name); i++){
		    if(key==options[i].key){
			iopt=i;
			break;
		    }
		}
	    }
	    if(iopt==-1){
		continue;/*we don't want this key. */
	    }
	    if((options[iopt].valtype)){//expects a value.
		value=start;
		while(value[0]=='\n' || isspace((int)value[0])){
		    value[0]=' ';
		    value++;
		}
	    }else{
		value=NULL;
	    }
	    int isfun=(options[iopt].isfun);
	    switch(options[iopt].type){
	    case 0:/*no result needed */
		break;
	    case M_INT:{
		if(options[iopt].valtype==2){//needs an array
		    if(isfun) error("Not implemented yet\n");
		    int val=strtol(value, &start, 10);
		    int **tmp=(int**)options[iopt].val;
		    int *nval=(int*)options[iopt].nval;
		    int i;
		    for(i=0; i<*nval; i++){
			if((*tmp)[i]==val) break;
		    }
		    if(i==*nval){
			(*nval)++;
			*tmp=myrealloc(*tmp, *nval,int);
			(*tmp)[(*nval)-1]=val;
		    } 
		}else{
		    int val=value?strtol(value, &start, 10):1;
		    if(isfun){/*Is function */
			void (*tmp)(int)=(void (*)(int))options[iopt].val;
			tmp(val);
		    }else{
			int *tmp=(int*)options[iopt].val;
			*tmp=val;
		    }
		}
	    }
		break;
	    case M_DBL:{
		if(options[iopt].valtype==2){//needs an array
		    if(isfun) error("Not implemented yet\n");
		    double val=strtod(value, &start);
		    double **tmp=(double**)options[iopt].val;
		    int *nval=(int*)options[iopt].nval;
		    (*nval)++;
		    *tmp=myrealloc(*tmp, *nval,double);
		    (*tmp)[(*nval)-1]=(int)val;
		}else{
		    double val=value?strtod(value, &start):1;
		    if(isfun){/*Is function */
			void (*tmp)(double)=(void (*)(double))options[iopt].val;
			tmp(val);
		    }else{
			double *tmp=(double*)options[iopt].val;
			*tmp=val;
		    }
		}
	    }
		break;
	    case M_STR:{
		char *val=value?cmd_string(value, &start):strdup("Unknown");
		if(isfun){
		    void (*tmp)(char*)=(void (*)(char*))options[iopt].val;
		    tmp(val);
		    free(val);
		}else{
		    char **tmp=(char**)options[iopt].val;
		    free(*tmp); *tmp=val;
		}
	    }
		break;
	    default:
		error("Unknown type");
	    }/*switch */
	    /*Empty the string that we already parsed. */
	    memset(start0, ' ',start-start0);
	}else if(start[0]=='='){/*equal sign found, key=value */
	    /*create a \n before the key. */
	    int skipspace=1;
	    for(char *start2=start-1; start2>=cmds; start2--){
  	        if(isspace((int)*start2) || *start2=='\n'){
		    if(!skipspace){
			*start2='\n';
			break;
		    }
		}else{
		    skipspace=0;
		}
	    }
	    start++;
	}else if(!mystrcmp(start, ".conf")){ /*.conf found. */
	    /*create a \n before the key. and \n after .conf */
	    for(char *start2=start-1; start2>=cmds; start2--){
	        if(isspace((int)*start2) || *start2=='\n'){
		    *start2='\n'; 
		    break;
		}
	    }
	    start+=5;
	    start[0]='\n';
	    start++;
	}else if(start[0]=='['){/*make sure we don't split brackets that are part of value. */
	    char *bend=strchr(start+1, ']');
	    char *bnextstart=strchr(start+1, '[');
	    if(bend && (!bnextstart || bend<bnextstart)){/*There is a closing bracket */
		for(; start<bend+1; start++){
		    if(start[0]=='\n') start[0]=' ';
		}
	    }else{
		error("Bracket is not closed\n");
		start++;
	    }
	}else if(start[0]=='\'' || start[0]=='"'){/*make sure we don't split strings that are part of value. */
	    char *quoteend=strchr(start, start[0]);
	    if(quoteend){
		start=quoteend+1;
	    }else{
		warning("Quote is not closed\n");
		start++;
	    }
	}else{
	    start++;
	}
    }
}
#include <semaphore.h>
/**
   Block signals in critical region.
*/
int sig_block(int block){
    sigset_t set;
    sigfillset(&set);
    if(block){
	return sigprocmask(SIG_BLOCK, &set, NULL);
    }else{
	return sigprocmask(SIG_UNBLOCK, &set, NULL);
    }
}
int sem_lock(const char *key, int lock){
    char fnsem[PATH_MAX];
    snprintf(fnsem, PATH_MAX, "%s_%s", key, USER);
    if(fnsem[0]!='/'){
	fnsem[0]='/';
    }
    sem_t *sem=sem_open(fnsem, lock==1?O_CREAT:0, 00700, 1);
    if(sem==SEM_FAILED){
	info("fnsem=%s\n", fnsem);
	perror("sem_open");
	return -1;
    }else{
	int value;
	sem_getvalue(sem, &value);
	if(lock){
	    info2("locking %p (value=%d) ... ", sem, value);
	    sem_wait(sem);
	}else{
	    info2("unlock %p (value=%d) ... ", sem, value);
	    sem_post(sem);
	}
	(void)sig_block(lock);
	info2("done\n");
	return 0;
    }
}
 
/**
   Set scheduling priorities for the process to enable real time behavior.
*/
void set_realtime(int icpu, int niceness){
    (void)icpu;
    (void)niceness;
    //Set CPU affinity.
    /*
      //Deprecated. Use external tools to do so, like openmp env's
#ifdef __linux__    
    if(icpu>0){
	cpu_set_t cpuset={{0}};
	CPU_SET(icpu, &cpuset);
	sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
    }
    //lock data in memory, avoid swapping.
    mlockall(MCL_FUTURE | MCL_CURRENT);
#endif
    */
    //fail stack
    struct rlimit rl;
    if(!getrlimit(RLIMIT_STACK, &rl)){
	const int NSTACK=rl.rlim_cur/2;
	char tmp[NSTACK];
	memset(tmp, 0, NSTACK);
    }
    //Set only if we are root.
    if(getuid()==0){
	info2("Set priority to -20\n");
	setpriority(PRIO_PROCESS, getpid(), -20);//is this necessary?
#ifdef __linux__
	struct sched_param param;
	sched_getparam(getpid(), &param);
	param.sched_priority=sched_get_priority_max(SCHED_FIFO)-1;
	sched_setscheduler(getpid(), SCHED_FIFO, &param);
#endif
    }else{
	warning("Please run program as setsid or as root to lift priority\n");
    }
}
quitfun_t quitfun=0;
void default_quitfun(const char *msg){
    fprintf(stderr, "%s", msg);
    sync();
    if(strncmp(msg, "ERROR", 5)){
	print_backtrace();
	sync();
    }
    exit(0);
}
static int (*signal_handler)(int)=0;
static volatile sig_atomic_t fatal_error_in_progress=0;
void default_signal_handler(int sig, siginfo_t *siginfo, void *unused){
    (void)unused;
    info2("default_signal_handler: %s (%d).\n", sys_siglist[sig], sig);sync();
    int cancel_action=0;
    struct sigaction act={{0}};
    act.sa_flags=0;
    act.sa_handler=SIG_DFL;
    sigaction(sig, &act, 0);
    /*prevent recursive call of handler*/
    if(sig==0){
	info2("Signal 0 caught. do nothing\n");
	sync();
	return;
    }
    if(fatal_error_in_progress){
	return;
    }
    extern int exit_fail;
    exit_fail=1;
    fatal_error_in_progress++;
    if(sig != SIGABRT && siginfo && siginfo->si_addr){
	info2("Memory location: %p\n", siginfo->si_addr);
    }
    if(sig==SIGBUS || sig==SIGILL || sig==SIGSEGV || sig==SIGABRT){
	print_backtrace();
    }
    if(signal_handler){
	info2("Signal %d caught. Call signal handler.\n", sig);
	if(signal_handler(sig)){
	    cancel_action=1; 
	}
    }else{
	info2("Signal %d caught without active handler.\n", sig);
    }
    sync();
    if(!cancel_action){//Propagate signal to default handler.
	act.sa_handler=SIG_DFL;
	sigaction(sig, &act, 0);
	raise(sig);
    }
}

/**
   Register signal handler
*/
void register_signal_handler(int (*func)(int)){
    struct sigaction act={{0}};
    act.sa_sigaction=default_signal_handler;
    act.sa_flags=SA_SIGINFO;
    sigaction(SIGBUS, &act, 0);
    sigaction(SIGILL, &act, 0);
    sigaction(SIGSEGV,&act, 0);
    sigaction(SIGINT, &act, 0);
    sigaction(SIGTERM,&act, 0);
    sigaction(SIGABRT,&act, 0);
    sigaction(SIGHUP, &act, 0);
    sigaction(SIGUSR1,&act, 0);
    sigaction(SIGQUIT,&act, 0);
    signal_handler=func;
}
/**
   Pause, waiting for user input
*/
void mypause(){
    info2("Press Any Key to Continue:"); 
    while(getchar()!=0x0a); 
    info2("continuing...\n"); 
}
#undef strdup
char* (*strdup0)(const char *)=strdup;

