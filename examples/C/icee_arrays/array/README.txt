
This is an example to demonstrate Adios ICEE method for WAN staging.

== ICEE Initialization ==

In short, ICEE is sending information in a push way. I.e., once a
reader registers his connection information (hostname and port) to a
writer, the writer pushes data to the reader.

In ICEE, an initial connection is established between a writer and a
reader by the following steps:

1. Writer listens on a port, say P1, waiting a client (or multiple) to
register. Writer's listening port can be specified by "cm_port"
option when calling adios_select_method (when using Adios no xml init)
or in an xml file.
and
The following options can be used for initializing writer:

    cm_port : writer's listening port
    max_client : number of clients to register (default is 1)

2. A reader needs to know the writer's hostname and port. Then, the
reader connects to the writer and registers his connection information
(reader's hostname and port). 

The following options can be specified on calling
adios_read_init_method:

    cm_host : reader's hostname
    cm_port : reader's port
    cm_remote_host : writer's hostname
    cm_remote_port : writer's port

A reader can register to multiple writers to receive data from
multiple sources. Then, use the following option:

    remote_list : specify multiple writer's location in the format:
    hostname1:port1,hostname2:port2,...

3. After successful registration, the writer will push data to the reader.



== ICEE Example ==

This program consists of two Adios executables; adios_write and
adios_read. adios_write will send data through TCP/IP to adios_read.

Case I: 1-to-1 connection (one writer and one reader)

We assume the following TCP/IP connection:

       adios_write            adios_read
node:    host1                  host2
port:    59997                  59999

Execution is as follows. First, run the writer:

$ adios_write -w ICEE -p 59997

Then, run the reader in another terminal:

$ adios_read -r ICEE -h host2 -p 59999 -s host1 -t 59997


Case II: N-to-1 connection (N writers and one reader)

We assume the following TCP/IP connection:

                  2 adios_write            adios_read
node:               host1                    host2
rank 0's port:      59997                    59999
rank 1's port:      59998                    

$ mpirun -n 2 adios_write -w ICEE -p 59997

Note: the rank 1's port will be computed by rank 0's port + 1. This is
not a ICEE function. It is programed in this example.

$ adios_read -r ICEE -h host2 -p 59999 -u "host1:59997,host1:59998"


Case III: 1-to-N connection (1 writer and N readers)

We assume the following TCP/IP connection:

                  adios_write            2 adios_read
node:               host1                    host2
rank 0's port:      59997                    59999
rank 1's port:                               60000


$ adios_write -w ICEE -p 59997 -m 2

$ mpirun -n 2 adios_read -r ICEE -h host2 -p 59999 -u "host1:59997"

Or, 
$ mpirun -n 2 adios_read -r ICEE -h host2 -p 59999 -s host1 -t 59997
