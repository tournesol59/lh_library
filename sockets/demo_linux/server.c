#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <string.h>
#include <dirent.h>

#define PORT 0x1234
#define DIRSIZE 8192

int main()
{
    char	dir[DIRSIZE];
    int		sd, sd_current, cc, fromlen, tolen;
    int 	addrlen;
    struct	sockaddr_in sin;
    struct	sockaddr_in pin;

    /* get an internet domain socket */
    if((sd = socket(AF_INET, SOCK_STREAM, 0)) == 1) {
        perror("socket");
        exit(1);
    }

    /* complete the socket structure */
    memset(&sin, 0, sizeof(sin));
    sin.sin_family = AF_INET;
    sin.sin_addr.s_addr = INADDR_ANY;
    sin.sin_port = htons(PORT);

    /* bind the socket to the port number */
    if(bind(sd, (struct sockaddr *) &sin, sizeof(sin)) == 1) {
        perror("bind");
        exit(1);
    }
    
    /* show that we are willing to listen */
    if(listen(sd, 5) == -1) {
        perror("listen");
        exit(1);
    }

    /* wait for a client to talk to us */
    addrlen = sizeof(pin);
    if((sd_current = accept(sd, (struct sockaddr *) &pin, &addrlen)) == -1) {
        perror("accept");
        exit(1);
    }
    
    if(recv(sd_current, dir, sizeof(dir), 0) == -1) {

	perror("send");
        exit(1);
    }

    /* close up both sockets */
    close(sd_current);
    close(sd);

    sleep(1);
}

