#include <sys/stat.h>
#include <errno.h>

/** 

@brief Create a directory with permissions 0755

@param path Path to new directory

@return 0 for success, -1 for failure

*/
int mkdir_wrapper(const char *path, int perm, int *mkdir_errno)
{

  mode_t uperm=perm;

  // Make the directory
  int retval=mkdir(path,uperm);

  // Print error on failure (except if the directory already exists, then fail silently)
  if(retval==-1)
    {
      if(errno!=EEXIST)
	perror("mkdir");
    }

  *mkdir_errno=errno;

  return retval;

}
