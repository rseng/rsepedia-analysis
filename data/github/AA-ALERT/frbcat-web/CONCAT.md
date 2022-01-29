[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1400274.svg)](https://doi.org/10.5281/zenodo.1400274)

# frbcat-web
This is the web frontend to the FRB Catalogue. This catalogue contains up to date information for the published population of Fast Radio Bursts (FRBs). The catalog connects to a postgres database, which is defined in the frbcatdb python package. For detailed information on the database format, see the frbcatdb package, or directly the [create table script](https://github.com/TRASAL/frbcatdb/blob/master/db/create_db_tables.sql).

### Configuration:
frbcat-web has two configuration files.

  * server/config.js

    This file contains the postgres connection string. Modify it to connect to your postgres installation and table. The connectionstring has the following format: 'postgres://someuser:somepassword@somehost:someport/somedatabase'
    
    In addition, the file contains a link to the CSV location at Zenodo, e.g. https://doi.org/10.5281/zenodo.xxxxxxxx
 
  * server.js
  
    This file contains the express web server configuration. Here you can define the port the server is running on.


For further details on how to tune the express server for performance in a production environment see the [expressjs homepage](https://expressjs.com/en/advanced/best-practice-performance.html). This link contains information on how to set the NODE_ENV to production for different systems and how to make sure the server gets automatically restarted if it crashes.

### Installation:
To install frbcat-web, please follow the following steps. First we need to install the dependencies of frbcat-web using the following commands:
```
  npm ci
```
Next, we build the package using:

```
  npm run webpack
```

### Running the express server with frbcat-web:
To start the express server with frbcat-web loaded use the following command:
```
  npm run start
```
