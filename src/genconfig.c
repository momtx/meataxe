#include <string.h>
#include <stdio.h>
#include <time.h>

#define MKSTRING2(x) #x
#define MKSTRING(x) MKSTRING2(x)

static int IsBigEndian()
{
    union {
        unsigned char c[sizeof(int)];
	unsigned u;
    } x;
    x.u = 1;
    return x.c[0] == 0;
}


static void printConfig()
{

    printf("#define MTX_CONFIG_BIG_ENDIAN %d\n",IsBigEndian());
    printf("#define MTX_CONFIG_LONG32 %d\n",sizeof(long) == 4);
    printf("#define MTX_CONFIG_LONG64 %d\n",sizeof(long) == 8);

    printf("#define MTX_CONFIG \"");
    switch (sizeof(long)) {
	case 8: printf("L64"); break;
	case 4: printf("L32"); break;
	default: printf("??"); break;
    }
    printf(IsBigEndian() ? " BE" : " LE");
    printf(" ZZZ=%d", ZZZ);
    printf("\"\n");

    time_t now = time(0);
    struct tm *tm = localtime(&now);
    printf("#define MTXBUILDTIME \"%04d-%02d-%02d/%0d:%02d:%02d\"\n", tm->tm_year + 1900,
	    tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec);
    printf("#define MTX_VERSION \"%s\"\n",MKSTRING(MTXVERSION) );
    printf("#define MTXVERSION \"%s\"\n",MKSTRING(MTXVERSION) );
}


int main()
{
    char line[1000];

    while (fgets(line, sizeof(line), stdin) != NULL) {
	if (strncmp(line, "@@insert_config_here", 20) == 0) {
	    printConfig();
	} else {
	    fputs(line,stdout);
	}
    }
}
