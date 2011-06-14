#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

//#include <arpa/inet.h>
//#include <sys/socket.h>
#include <netdb.h>



int main(void) 
{

    int sockfd;  
    struct addrinfo hints, *servinfo, *p;
    int rv;
    char host[NI_MAXHOST];

    memset(&hints, 0, sizeof hints);
    hints.ai_family = AF_UNSPEC; // use AF_INET6 to force IPv6
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE; // use my IP address

    if ((rv = getaddrinfo(NULL, "34909", &hints, &servinfo)) != 0) {
        fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(rv));
        exit(1);
    }

    // loop through all the results and bind to the first we can
    for(p = servinfo; p != NULL; p = p->ai_next) {
        /*
        if ((sockfd = socket(p->ai_family, p->ai_socktype,
                        p->ai_protocol)) == -1) {
            perror("socket");
            continue;
        }

        if (bind(sockfd, p->ai_addr, p->ai_addrlen) == -1) {
            close(sockfd);
            perror("bind");
            continue;
        }

        break; // if we get here, we must have connected successfully
        */
        int s = getnameinfo(p->ai_addr, sizeof(struct sockaddr_in),
                host, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);
        if (s != 0) {
            printf("getnameinfo() failed: %s\n", gai_strerror(s));
            exit(EXIT_FAILURE);
        }
        //printf("<Interface>: %s \t <Address> %s\n", ifa->ifa_name, host);


        printf ("Name=%s, host=%s, flags=%d, family=%d, socktype=%d, protocol=%d\n",
               p->ai_canonname, host, p->ai_flags, p->ai_family, p->ai_socktype, p->ai_protocol);
    }

    /*
    if (p == NULL) {
        // looped off the end of the list with no successful bind
        fprintf(stderr, "failed to bind socket\n");
        exit(2);
    }
    */

    freeaddrinfo(servinfo); // all done with this structure
}
