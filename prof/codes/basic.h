#ifndef _BASIC_
#define _BASIC_

void my_error(char *error_text) ;
void indexxd(int n, double arr[], int indx[]);
int mkdirectory(char *dir);
int test_file(char *dir);
int check_parent(char *dir);
int nofile(FILE *fp, char *filename);
int mk_savedir(char *pre);

#endif

