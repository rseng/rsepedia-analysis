# Elixir Daisy
[![Build Status](https://travis-ci.com/elixir-luxembourg/daisy.svg?branch=develop)](https://travis-ci.com/elixir-luxembourg/daisy)

Data Information System (DAISY) is a data bookkeeping application designed to help Biomedical Research institutions with their GDPR compliance.

For more information, please refer to the official [Daisy documentation](https://elixir.pages.uni.lu/daisy-doc/).

## Demo deployment
You are encouraged to try Daisy for yourself using our [DEMO deployment](https://daisy-demo.elixir-luxembourg.org/).

## Deployment using Docker

### Requirements

* docker: https://docs.docker.com/install/

### Installation

1. Get the source code
    
    ```bash
    git clone git@github.com:elixir-luxembourg/daisy.git
    cd daisy
    ```
1. Create your settings file
    
	```bash
	cp elixir_daisy/settings_local.template.py elixir_daisy/settings_local.py
	```
    Optional: edit the file elixir_daisy/settings_local.py to adapt to your environment.

1. Build daisy docker image  
    ```bash
    docker-compose up --build
    ```
    Wait for the build to finish and keep the process running
1. Open a new shell and go to daisy folder

1. Build the database
    
    ```bash
    docker-compose exec web python manage.py migrate
    ```
1. Build the solr schema

    ```bash
    docker-compose exec web python manage.py build_solr_schema -c /solr/daisy/conf -r daisy -u default
    ```

1. Compile and deploy static files
    
    ```bash
    docker-compose exec web python manage.py collectstatic
    ```
1. Create initial data in the database
    
    ```bash
    docker-compose exec web bash -c "cd core/fixtures/ && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/edda.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hpo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hdo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hgnc.json"
    docker-compose exec web python manage.py load_initial_data
    ```
   Initial data includes, for instance, controlled vocabularies terms and initial list of institutions and cohorts.  
   **This step can take several minutes to complete**
    
1. Load demo data
    
    ```bash
    docker-compose exec web python manage.py load_demo_data
    ```
    This will create mock datasets, projects and create an demo admin account.

1. Optional - import users from an active directory instance

    ```bash
    docker-compose exec web python manage.py import_users
    ```
    
1.  Build the search index

    ```bash
    docker-compose exec web python manage.py rebuild_index -u default
    ```    

1. Browse to https://localhost  
    a demo admin account is available:
    
    ```
        username: admin
        password: demo
    ```

### Operation manual


#### Importing 

In addition to loading of initial data, DAISY database can be populated by importing Project, Dataset and Partners records from JSON files using commands `import_projects`, `import_datasets` and `import_partners` respectively.
 The commands for import are accepting one JSON file (flag `-f`): </br>

```bash
docker-compose exec web python manage.py <COMMAND> -f ${PATH_TO_JSON_FILE}
```
where ${PATH_TO_JSON_FILE} is the path to a json file containing the records definitions.
See file daisy/data/demo/projects.json as an example.
 
Alternatively, you can specify directory containing multiple JSON files to be imported with `-d` flag:
```bash
docker-compose exec web python manage.py <COMMAND> -d ${PATH_TO_DIR}
```

#### Exporting  

Information in the DAISY database can be exported to JSON files. The command for export are given below:</br>

```bash
docker-compose exec web python manage.py export_partners -f ${JSON_FILE}
```
where ${JSON_FILE} is the path to a json file that will be produced.  In addition to ````export_partners````, you can run ````export_projects```` and ````export_datasets```` in the same way.

### Upgrade to last Daisy version

1. Create a database backup.

	```bash
	docker-compose exec db pg_dump daisy --port=5432 --username=daisy --no-password --clean > backup_`date +%y-%m-%d`.sql
	```
	
1. Make sure docker containers are stopped.

	```bash
	docker-compose stop
	```

3. Get last Daisy release.

	```bash
	git checkout master
	git pull
	```

1. Rebuild and start the docker containers.

	```bash
	docker-compose up --build
	```
	Open a new terminal window to execute the following commands.

1. Update the database schema.

	```bash
	docker-compose exec web python manage.py migrate
	```

1. Update the solr schema.

	```bash
	docker-compose exec web python manage.py build_solr_schema -c /solr/daisy/conf -r daisy -u default
	```

1. Collect static files.
 
	```bash
	docker-compose exec web python manage.py collectstatic
	```
	
1. Reload initial data (optional).


    **IMPORTANT NOTE:** The initial data package provides some default values for various lookup lists e.g. data sensitivity classes, document or data types.  If, while using DAISY, you have customized these default lists, please keep in mind that running the ```load_initial_data``` command during update will re-introduce those default values. If this is not desired, then please skip the reloading of initial data step during your update. You manage lookup lists through the application interface.<br/><br/>

 
    ```bash
    docker-compose exec web bash -c "cd core/fixtures/ && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/edda.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hpo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hdo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hgnc.json"
    docker-compose exec web python manage.py load_initial_data
    ```
	
	
    **IMPORTANT NOTE:** This step can take several minutes to complete. 
	
1. Rebuild the search index.
    
    ```bash
    docker-compose exec web python manage.py rebuild_index -u default
    ```	
1. Reimport the users (optional).
	    
    If LDAP was used during initial setup to import users, they have to be imported again:
    
    ```bash
    docker-compose exec web python manage.py import_users
    ```
    
## Deployment without Docker - CentOS


See [DEPLOYMENT](DEPLOYMENT.md).


## Development

To be completed.

### Import users from active directory
```bash
./manage.py import_users
```

### Import projects, datasets or partners from external system
Single file mode:
```bash
./manage.py import_projects -f path/to/json_file.json
```

Batch mode:
```bash
./manage.py import_projects -d path/to/dir/with/json/files/
```

Available commands: `import_projects`, `import_datasets`, `import_partners`.

In case of problems, add `--verbose` flag to the command, and take a look inside `./log/daisy.log`. 

### Install js and css dependencies

```bash
cd web/static/vendor/
npm ci
```

### Compile daisy.scss
```bash
cd web/static/vendor
npm run-script build
```

### Run the built-in web server (for development)

```bash
./manage.py runserver
```

### Run the tests

The following command will install the test dependencies and execute the tests:

```bash
python setup.py pytest
```

If tests dependencies are already installed, one can also run the tests just by executing:

```bash
pytest
```

## Administration

To get access to the admin page, you must log in with a superuser account.  
On the `Users` section, you can give any user a `staff` status and he will be able to access any project/datasets.


## `settings.py` and `local_settings.py` reference

### Display
| Key | Description | Expected values | Example value |
|---|---|---|---|
| `COMPANY` | A name that is used to generate verbose names of some models | str | `'LCSB'` |
| `DEMO_MODE` | A flag which makes a simple banneer about demo mode appear in About page | bool | `False` | 
| `INSTANCE_LABEL` | A name that is used in navbar header to help differentiate different deployments | str | `'Staging test VM'` | 
| `INSTANCE_PRIMARY_COLOR` | A color that will be navbar header's background | str of a color | `'#076505'` | 
| `LOGIN_USERNAME_PLACEHOLDER` | A helpful placeholder in login form for logins | str | `'@uni.lu'` | 
| `LOGIN_PASSWORD_PLACEHOLDER` | A helpful placeholder in login form for passwords | str | `'Hint: use your AD password'` | 

### Integration with other tools
#### ID Service
| Key | Description | Expected values | Example value |
|---|---|---|---|
| `IDSERVICE_FUNCTION` | Path to a function (`lambda: str`) that generates IDs for entities which are published | str | `'web.views.utils.generate_elu_accession'` |
| `IDSERVICE_ENDPOINT` | In case LCSB's idservice function is being used, the setting contains the IDservice's URI | str | `'https://192.168.1.101/v1/api/` |

#### REMS
| Key | Description | Expected values | Example value |
|---|---|---|---|
| `REMS_INTEGRATION_ENABLED` | A feature flag for REMS integration. In practice, there's a dedicated endpoint which processes the information from REMS about dataset entitlements | str | `True` |
| `REMS_MATCH_USERS_BY` | A method of how the information from REMS is mapped to the Users - via emails (`email`) or OIDC ID (`id`), or both (`auto`) | str['auto','email','id'] | `False` |
| `REMS_SKIP_IP_CHECK` | If set to `True`, there will be no IP checking if the request comes from trusted REMS instance. | bool | `False` |
| `REMS_ALLOWED_IP_ADDRESSES` | A list of IP addresses that should be considered trusted REMS instances. Beware of configuration difficulties when using reverse proxies. The check can be skipped with `REMS_SKIP_IP_CHECK` | dict[str] | `['127.0.0.1', '192.168.1.101']` |

#### Keycloak
| Key | Description | Expected values | Example value |
|---|---|---|---|
| `KEYCLOAK_INTEGRATION` | A feature flag for importing user information from Keycloak (OIDC IDs) | bool | `True` |
| `KEYCLOAK_URL` | If keycloak integration flag is set, this setting contains URL to the Keycloak instance | str | `'https://keycloak.lcsb.uni.lu/auth/'` |
| `KEYCLOAK_REALM` | If keycloak integration flag is set, this setting contains the realm's name in your Keycloak instance | str | `'master'` |
| `KEYCLOAK_USER` | If keycloak integration flag is set, this setting contains the username to access Keycloak | str | `'username'` |
| `KEYCLOAK_PASS` | If keycloak integration flag is set, this setting contains the password to access Keycloak | str | `'secure123'` |

### Others
| Key | Description | Expected values | Example value |
|---|---|---|---|
| `SERVER_SCHEME` | A URL's scheme to access your DAISY instance (http or https) | str | `'https'` |
| `SERVER_URL` | A URL to access your DAISY instance (without the scheme) | str | `'example.com'` |
| `GLOBAL_API_KEY` | An API key that is not connected with any user. Disabled if set to `None` | optional[str] | `'in-practice-you-dont-want-to-use-it-unless-debugging'` |# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at info@elixir-luxembourg.org. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Installation

## Base

```bash
sudo yum update
sudo yum install python36-devel openldap-devel nginx
sudo yum group install "Development Tools"

wget https://bootstrap.pypa.io/get-pip.py
sudo python3.6 get-pip.py
```



# User and Application Source Code

We install the application under a dedicated `daisy` user.

```bash
sudo useradd daisy
sudo usermod -a -G users daisy
sudo su - daisy
mkdir config log
git clone https://github.com/elixir-luxembourg/daisy.git
exit
sudo /usr/local/bin/pip install -e /home/daisy/daisy
sudo /usr/local/bin/pip install gunicorn
```

## NPM and Node.js

```bash
curl -sL https://rpm.nodesource.com/setup_10.x | sudo bash -
sudo yum install nodejs
```

Then you need to compile the static files.


```bash
sudo su - daisy
cd /home/daisy/daisy/web/static/vendor/
npm ci
exit
```


## Solr

```bash
sudo useradd solr
wget https://archive.apache.org/dist/lucene/solr/7.7.1/solr-7.7.1.tgz
tar -xf solr-7.7.1.tgz solr-7.7.1/bin/install_solr_service.sh
sudo yum install lsof java-1.8.0-openjdk
sudo solr-7.7.1/bin/install_solr_service.sh solr-7.7.1.tgz
sudo su - solr
/opt/solr-7.7.1/bin/solr create_core -c daisy
cd /var/solr/data/daisy/conf
wget "https://raw.githubusercontent.com/apache/lucene-solr/master/solr/example/files/conf/currency.xml"  
wget "https://raw.githubusercontent.com/apache/lucene-solr/master/solr/example/files/conf/elevate.xml"  
/opt/solr-7.7.1/bin/solr stop  
exit
```
It is possible that by this time solr-7.7.1 is not anymore proposed for download on solr mirrors.
In this case check for last solr version available and adapt the instructions above accordingly.
You need configure the solr core 'daisy'. To do so you need to create 'schema.xml' and 'solrconfig.xml' files under 
'/var/solr/data/daisy/conf'. 
```bash
sudo cp /home/daisy/daisy/docker/solr/schema.xml /var/solr/data/daisy/conf/
sudo cp /home/daisy/daisy/docker/solr/solrconfig.xml /var/solr/data/daisy/conf/
```

Grant ownership and change privileges of `/var/solr` folder
```
sudo chown -R solr:users /var/solr
sudo chmod -R 775 /var/solr
```

<span style="color:red;">Review the 'schema.xml' file you just copied. Ensure that all file references inside it (e.g. stopwords.txt) actually exist in the path specified.</span>

<p style="color:red;">By default, the Solr instance listens on port 8983 on all interfaces.
Solr has no authentication system. It is crucial to secure it by either blocking external accesses to the Solr port or by changing it's configuration to listen only on localhost (see https://stackoverflow.com/a/1955591)</p>


You can restart solr and check that it is working with the following commands

```bash
sudo systemctl enable solr
sudo systemctl restart solr
```

## Gunicorn

1) Create the file ```/etc/systemd/system/gunicorn.service``` as the _root_ user or with _sudo_ and with the following content:


```
[Unit]
Description=gunicorn daemon
After=network.target

[Service]
PIDFile=/run/gunicorn/pid
User=daisy
Group=daisy
WorkingDirectory=/home/daisy/daisy
ExecStart=/usr/local/bin/gunicorn --limit-request-line 0 --access-logfile /home/daisy/log/gunicorn_access.log --error-logfile /home/daisy/log/gunicorn_error.log --log-level debug --workers 2 --bind unix:/home/daisy/daisy/daisy.sock elixir_daisy.wsgi
ExecReload=/bin/kill -s HUP $MAINPID
ExecStop=/bin/kill -s TERM  $MAINPID

[Install]
WantedBy=multi-user.target
```

## Rabbitmq

```bash
sudo yum install rabbitmq-server  
sudo systemctl start rabbitmq-server  
sudo systemctl enable gunicorn  
```

## Celery

We use systemd to create two services, celery_worker to run the worker (notifications, indexation, etc) and celery_beat to run the scheduled tasks.

1) Celery worker

As daisy user, create the file /home/daisy/config/celery.conf with the following content:

```
# Name of nodes to start
# here we have a single node
CELERYD_NODES="daisy_worker"
# or we could have three nodes:
#CELERYD_NODES="w1 w2 w3"

# Absolute or relative path to the 'celery' command:
CELERY_BIN="/usr/local/bin/celery"
#CELERY_BIN="/virtualenvs/def/bin/celery"

# App instance to use
# comment out this line if you don't use an app
CELERY_APP="elixir_daisy.celery_app"
# or fully qualified:
#CELERY_APP="proj.tasks:app"

# How to call manage.py
CELERYD_MULTI="multi"

# Extra command-line arguments to the worker
CELERYD_OPTS="--concurrency=1"

# - %n will be replaced with the first part of the nodename.
# - %I will be replaced with the current child process index
#   and is important when using the prefork pool to avoid race conditions.
CELERYD_PID_FILE="/var/run/celery/%n.pid"
CELERYD_LOG_FILE="/home/daisy/log/celery/%n%I.log"
CELERYD_LOG_LEVEL="DEBUG"
```

Create the  folders '/var/run/celery/' as _root_ or with _sudo_ and the folder '/home/daisy/log/celery' as _daisy_ must be created. 
Create also the service config file '/etc/systemd/system/celery_worker.service' as _root_ or with _sudo_ and with the following content:

```
[Unit]
Description=Celery Worker
After=network.target

[Service]
Type=forking
User=daisy
Group=daisy
EnvironmentFile=/home/daisy/config/celery.conf
WorkingDirectory=/home/daisy/daisy
ExecStart=/bin/sh -c '${CELERY_BIN} multi start ${CELERYD_NODES} \
  -A ${CELERY_APP} --pidfile=${CELERYD_PID_FILE} \
  --logfile=${CELERYD_LOG_FILE} --loglevel=${CELERYD_LOG_LEVEL} ${CELERYD_OPTS}'
ExecStop=/bin/sh -c '${CELERY_BIN} multi stopwait ${CELERYD_NODES} \
  --pidfile=${CELERYD_PID_FILE}'
ExecReload=/bin/sh -c '${CELERY_BIN} multi restart ${CELERYD_NODES} \
  -A ${CELERY_APP} --pidfile=${CELERYD_PID_FILE} \
  --logfile=${CELERYD_LOG_FILE} --loglevel=${CELERYD_LOG_LEVEL} ${CELERYD_OPTS}'

[Install]
WantedBy=multi-user.target
```

Then do the following:

```bash
sudo systemctl enable celery_worker  
sudo systemctl start celery_worker  
```

2) Celery beat

Create the file /home/daisy/config/celerybeat.conf with the following content:

```
# Absolute or relative path to the 'celery' command:
CELERY_BIN="/usr/local/bin/celery"

# App instance to use
# comment out this line if you don't use an app
CELERY_APP="elixir_daisy.celery_app"
# or fully qualified:
#CELERY_APP="proj.tasks:app"

# Extra command-line arguments to the worker
CELERYBEAT_OPTS="--scheduler django_celery_beat.schedulers:DatabaseScheduler"

# - %n will be replaced with the first part of the nodename.
# - %I will be replaced with the current child process index
#   and is important when using the prefork pool to avoid race conditions.
CELERYBEAT_PID_FILE="/var/run/celerybeat/%n.pid"
CELERYBEAT_LOG_FILE="/home/daisy/log/celerybeat/%n%I.log"
CELERYBEAT_LOG_LEVEL="INFO"
```

Create the service file /etc/systemd/system/celery_beat.service:

```
[Unit]
Description=Celery Beat Service
After=network.target

[Service]
User=daisy
Group=daisy
EnvironmentFile=/home/daisy/config/celerybeat.conf
WorkingDirectory=/home/daisy/daisy
ExecStart=/bin/sh -c '${CELERY_BIN} beat -A ${CELERY_APP} --pidfile=${CELERYBEAT_PID_FILE} \
  --logfile=${CELERYBEAT_LOG_FILE} --loglevel=${CELERYBEAT_LOG_LEVEL} ${CELERYBEAT_OPTS}'
ExecStop=/bin/kill -s TERM $MAINPID

[Install]
WantedBy=multi-user.target 
```

Then do the following:

```bash
sudo systemctl enable celery_beat 
sudo systemctl start celery_beat  
```

## PostgreSQL

### Install database server

```bash
sudo yum install https://download.postgresql.org/pub/repos/yum/10/redhat/rhel-7-x86_64/pgdg-centos10-10-2.noarch.rpm
sudo yum install postgresql10
sudo yum install postgresql10-server
sudo /usr/pgsql-10/bin/postgresql-10-setup initdb
sudo systemctl enable postgresql-10
sudo systemctl start postgresql-10
```

In case the installation fails, follow steps in the [official documentation](https://www.postgresql.org/download/linux/redhat/) for installation of Postgresql 10 on your platform.

### Create database and roles


```bash
sudo su - postgres
vi ./10/data/pg_hba.conf
```
Change METHOD ident of IPv4 and IPv6 to md5 and add rule for daisy and postgres users.
We recommend to only allow local connection from the daisy user to the daisy database.  
Example:  
```
# TYPE  DATABASE        USER            ADDRESS                 METHOD
local   all             postgres                                peer
# "local" is for Unix domain socket connections only
local   daisy           daisy                                   ident
local   postgres        postgres                                ident
# IPv4 local connections:
host    all             all             127.0.0.1/32            md5
# IPv6 local connections:
host    all             all             ::1/128                 md5
```

Create daisy user and database:
```
createuser daisy
createdb daisy
psql
postgres=# alter user daisy with encrypted password 'daisy';
postgres=# grant all privileges on database daisy to daisy ;
postgres=# \q
exit
```
<span style="color:red;">You can replace password `daisy` by a password of your choice.</span>

Restart PostgreSQL:

```bash
sudo systemctl restart postgresql-10
```

# Application Config File
Create a local configuration file for the application.

```bash
sudo su - daisy
cp /home/daisy/daisy/elixir_daisy/settings_local.template.py   /home/daisy/daisy/elixir_daisy/settings_local.py
vi /home/daisy/daisy/elixir_daisy/settings_local.py
```

<span style="color:red;">Change SECRET_KEY variable:</span>

```
# SECURITY WARNING: change the secret key used in production and keep it secret !
SECRET_KEY='<your-new-secret-key>'
```

Put in the following database configuration to the 'settings_local.py' file.

```
#......
#......
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': 'daisy',
        'USER': 'daisy',
        'PASSWORD': 'daisy',
        'HOST': 'localhost',
        'PORT': 5432
    }
}
#......
#......
```

Put in the following haystack configuration to the 'settings_local.py' file.

```
#......
#......
HAYSTACK_CONNECTIONS = {
        'default': {
            'ENGINE': 'haystack.backends.solr_backend.SolrEngine',
            'URL': 'http://127.0.0.1:8983/solr/daisy',
            'ADMIN_URL': 'http://127.0.0.1:8983/solr/admin/cores',
        },
}
HAYSTACK_SIGNAL_PROCESSOR = 'celery_haystack.signals.CelerySignalProcessor'
#......
#......
```
Add the following entries:
```
STATIC_ROOT = "/home/daisy/static/"
ALLOWED_HOSTS = ['10.x.x.x','daisy.com']
DEBUG = False
SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
```

Please note that the IP and DNS record should be CHANGED to denote your server.

Replace the company name 'LCSB' with your institution name.
We suggest that you use a not very long name here e.g. the acronym of your institution.

If needed, configure the active directory parameters to allow LDAP authentication and user imports. 
Exit the daisy user.
```bash
exit
```
# Web server

1) Install nginx

    ```bash
    sudo yum install epel-release
    sudo yum install nginx
    sudo systemctl enable nginx
    sudo systemctl start nginx
    ```
    
2) As _root_ or with _sudo_ create the file ```/etc/nginx/conf.d/ssl.conf``` with the following content:
	   
    ```bash
    proxy_connect_timeout       600;
    proxy_send_timeout          600;
    proxy_read_timeout          600;
    send_timeout                600;
    
    server {
        server_name daisy.com;
    
        location /static {
            alias /home/daisy/static;
            autoindex on;
        }
        
        location / {
            proxy_set_header Host $http_host;
            proxy_set_header X-Real-IP $remote_addr;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header X-Forwarded-Proto $scheme;
            proxy_pass http://unix:/home/daisy/daisy/daisy.sock;
        }
        listen 443 http2 ssl;
        listen [::]:443 http2 ssl;
        ssl on;
        ssl_certificate /etc/ssl/certs/daisy.com.crt;
        ssl_certificate_key /etc/ssl/private/daisy.com.key;
    }
	```    
    
	Changing daisy.com to your particular case.  
    
3) To have a redirect from http to https, as _root_ or with _sudo_ create the file ```/etc/nginx/conf.d/daisy.conf``` with the following content:

    ```
    server {
      listen 80;
      server_name daisy.com;
      return 301 https://daisy.com$request_uri;
    }
    ```
    Changing daisy.com to your particular case.
    
4) Create self-signed certificates if they already don't exist.

    ```bash
    openssl req -x509 -newkey rsa:4096 -nodes -out daisy.com.crt -keyout daisy.com.key -days 365
    ```
    Changing daisy.com to your particular case.
    Certificates should be put in the folder specified in ```/etc/nginx/conf.d/daisy.conf```
    ```bash
    sudo cp daisy.com.crt /etc/ssl/certs/
    sudo mkdir /etc/ssl/private/
    sudo cp daisy.com.key /etc/ssl/private/
 
    ```
5) Edit the file /etc/nginx/nginx.conf:

    Comment out the block server {} in /etc/nginx/nginx.conf
    Change the user running nginx from nginx to daisy
    
6) Grant access on `/var/lib/nginx` to **daisy** user:
   ```
   sudo chown -R daisy:daisy /var/lib/nginx
   ```
   
7) Restart nginx

    ```bash
    sudo systemctl restart nginx
    ```

# Initialization

Once everything is set up, the definitions and lookup values need to be inserted into the database.   
To do this run the following.

```bash
sudo su - daisy
cd /home/daisy/daisy
python3.6 manage.py collectstatic 
python3.6 manage.py migrate 
python3.6 manage.py build_solr_schema -c /var/solr/data/daisy/conf -r daisy  
cd /home/daisy/daisy/core/fixtures/
wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/edda.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hpo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hdo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hgnc.json
cd /home/daisy/daisy
python3.6 manage.py load_initial_data
```
The load_initial_data command needs several minutes to complete.
DAISY has a demo data loader. With example records of Projects Datasets and Users. If you want to deploy DAISY  demo data, then do 

```bash
python3.6 manage.py load_demo_data
```

The above command will create an 'admin' and other users such as 'alice.white', 'john.doe' 'jane.doe'. The password for all  is 'demo'.


If you do not want to load the demo data and work with your own definitions, then you'd still need to create super user for the application, with which you can logon and create other users as well as records. To create a super user, do the following and respond to the questions. 

```bash
python3.6 manage.py createsuperuser
```

Trigger a reindex with:

```bash
python3.6 manage.py rebuild_index
```

# Validate the installation

Check the the installation was successful by accessing the URL `https://${IP_OF_THE_SERVER}` with a web browser.
You should be able to login with `admin/demo` if the `load_demo_data` command was used or with your own admin account if the `createsuperuser` command was used.
It should be possible to create datasets and projects.

In addition when the DAISY is updated or configurations are changed (including the configuration files such as ```settings_local.py```) is modified, gunicorn must be restarted to load the new code/configuration, to do so run:

```bash
sudo systemctl restart gunicorn
sudo systemctl restart celery_worker
sudo systemctl restart celery_beat
```

# Setting up reminders 

DAISY can generate reminders on approaching deadlines (e.g. data storage end date or document expiry). To enable this feature, do the following:

 1) Login to DAISY as a super user. e.g. `admin` user in the demo application
 
 2) Go to https://${IP_OF_THE_SERVER}/admin
 
 3) From the 'Site administration' list select 'Periodic tasks' under 'PERIODIC TASKS' heading.
 
 4) Clicking the 'ADD PERIODIC TASK' button, then:
    4.1) Give your task a name, 
    4.2) From the 'Task(registered)' list select `notification.tasks.document_expiry_notifications`,
    4.3) From the 'Interval' list select `every day`. If this interval does not exist, you may create it by clicking the (+) button next to the select,.
    4.4) Select a start date and time, e.g. today and now,
    4.5) Click 'SAVE'.
    
 5) You may repeat the steps in (4) to create a daily periodic task also for `notification.tasks.data_storage_expiry_notifications`,
 
# Updating DAISY


If you want to move to the newest release of DAISY, do the following.  

1) Stop services, create a database and application backup.

As root user:

```bash
systemctl stop gunicorn
systemctl stop celery_worker
systemctl stop celery_beat 
su -c 'PGPASSWORD="<PASSWORD_OF_POSTGRES_USER>" pg_dump daisy --port=5432 --username=daisy --clean > daisy_dump.sql' - daisy 
tar -cvf /tmp/daisy.tar /home/daisy 
```

Once you have have created the tar ball of the application directory and the postgres dump, then you may proceed to update.

2) Get the latest Daisy release.

As daisy user:

```bash
cd /home/daisy/daisy
git checkout -- web/static/vendor/package-lock.json
git checkout master
git pull


cd /home/daisy/daisy/web/static/vendor/
npm ci
```
As root user:

```bash
/usr/local/bin/pip install -e /home/daisy/daisy --upgrade
```

3) Update the database and solr schemas, collect static files.

As daisy user:

```bash
cd /home/daisy/daisy
python3.6 manage.py migrate && python3.6 manage.py build_solr_schema -c /var/solr/data/daisy/conf/ -r daisy && yes | python3.6 manage.py clear_index && yes "yes" | python3.6 manage.py collectstatic;
```


4) Reload initial data (optional). 


**IMPORTANT NOTE:** The initial data package provides some default values for various lookup lists e.g. data sensitivity classes, document or data types.  If, while using DAISY, you have customized these default lists, please keep in mind that running the ``load_initial_data`` command
during update will re-introduce those default values. If this is not desired, then please skip the reloading of initial data step during your update. You manage lookup lists through the application interface.<br/><br/>
        
As daisy user:
```bash
cd /home/daisy/daisy/core/fixtures/
wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/edda.json -O edda.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hpo.json -O hpo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hdo.json -O hdo.json && wget https://git-r3lab.uni.lu/pinar.alper/metadata-tools/raw/master/metadata_tools/resources/hgnc.json -O hgnc.json

cd /home/daisy/daisy
python3.6 manage.py load_initial_data
```     
  
**IMPORTANT NOTE:** This step can take several minutes to complete. 


5) Reimport the users (optional).

If LDAP was used to import users, they have to be imported again.
As daisy user:
```bash
python3.6 manage.py import_users
```

6) Rebuild Solr search index.

As daisy user:
 ```bash
 cd /home/daisy/daisy
 python3.6 manage.py rebuild_index
```

7) Restart services.

As root user:

```bash
systemctl start gunicorn
systemctl start celery_worker
systemctl start celery_beat 
```
# Restoring backup of Daisy
First, make sure you have successfully backed up your Daisy deployment - see first section of chapter Updating Daisy.
Your backup .tar file should contain both the dump of Postgresql database and everything from `/home/daisy` directory.

As root user, stop services:
```bash
systemctl stop gunicorn
systemctl stop celery_worker
systemctl stop celery_beat
```

Wipe out broken/unwanted version of Daisy by deleting all files in daisy user home directory and dropping the database:  
**IMPORTANT NOTE**: Be sure that your backup .tar file is stored somewhere else!
```
rm -rf /home/daisy/*
su -c 'dropdb daisy' - postgres
```

Restore files from tar ball:
```
su -c 'tar -xvf <PATH-TO-BACKUP-FOLDER>/daisy.tar --directory /' - daisy
```

Following steps assume that the Postgresql10 is installed, pg_hba.conf file is updated and database user *daisy* exists (please see the postgresql deployment instructions for more information).
Create the database and grant privileges:
```
su - postgres
createdb daisy
psql -d daisy -p 5432 -c "grant all privileges on database daisy to daisy"
exit
```

Restore the database as daisy user: 
```
su -c 'psql -d daisy -U daisy -p 5432 < /home/daisy/daisy_dump.sql' - daisy
```

Start services: 
```
systemctl start gunicorn
systemctl start celery_worker
systemctl start celery_beat 
```
# bootstrap-material-datetimepicker
Material DateTimePicker

Originaly designed for Bootstrap Material, the V2.0 is now completely standalone and responsive.

### Updates

| Date				| Author			| Description											 |
| ----------------- | ----------------- | ------------------------------------------------------ |
| 2016-04-08		| donovansolms		| Disable specific days (#60 and #97)				 	 |
| 2016-04-08		| T00rk				| Fixed #85	 								 	 		 |
| 2016-04-08		| FoxyCorndog		| Fix PM overwrite bug	 					 	 		 |
| 2016-02-17		| T00rk				| Changed Clock to SVG	 					 	 		 |
| 2016-02-04		| T00rk				| Added a "Now" button (#38)	 					 	 |
| 2016-01-30		| T00rk				| Switch view on click (#39, #47)	 					 |
| 2016-01-29		| T00rk				| Added "clear button" (#48)		 					 |
| 2016-01-29		| T00rk				| Replace rem by em (#26)			 					 |
| 2016-01-29		| T00rk				| Display 24H clock (#54)			 					 |
| 2016-01-29		| T00rk				| Close on "ESC" (#52)			 					 	 |
| 2015-10-19		| drblue 			| Fixed erroneous package.json-file 					 |
| 2015-10-19		| Perdona			| Fix auto resize when month has 6 weeks				 |
| 2015-07-01		| T00rk 			| Redesigned element without using modal				 |
| 2015-06-16		| T00rk 			| Use Timepicker alone / Display short time (12 hours)	 |
| 2015-06-13		| T00rk 			| Fixed issue with HTML value tag 						 |
| 2015-05-25		| T00rk 			| Changed repo name to bootstrap-material-datetimepicker * |
| 2015-05-12		| T00rk				| Added parameters for button text						 |
| 2015-05-05		| Sovanna			| FIX undefined _minDate in isBeforeMaxDate func		 |
| 2015-04-10		| T00rk				| Little change in clock design							 |
| 2015-04-10		| Peterzen			| Added bower and requirejs support						 |
| 2015-04-08		| T00rk				| Fixed problem on locale switch						 |
| 2015-03-04		| T00rk				| Added Time picker										 |
(\*) File names have been changed

bootstrap-material-datepicker.js => bootstrap-material-date**time**picker.js

bootstrap-material-datepicker.css => bootstrap-material-date**time**picker.css

### Prerequisites

jquery [http://jquery.com/download/](http://jquery.com/download/)

momentjs [http://momentjs.com/](http://momentjs.com/)

Google Material Icon Font `<link href="https://fonts.googleapis.com/icon?family=Material+Icons" rel="stylesheet">`


### Live Example

[Live example](http://t00rk.github.io/bootstrap-material-datetimepicker/)

### Usage

	$('input').bootstrapMaterialDatePicker();

### bower

	bower install bootstrap-material-datetimepicker

### Parameters

| Name				| Type							| Description									|
| ----------------- | ----------------------------- | --------------------------------------------- |
| **format**		| String						| MomentJS Format								|
| **shortTime**		| Boolean						| true => Display 12 hours AM|PM 				|
| **minDate**		| (String\|Date\|Moment)		| Minimum selectable date						|
| **maxDate**		| (String\|Date\|Moment)		| Maximum selectable date						|
| **currentDate**	| (String\|Date\|Moment)		| Initial Date									|
| **year**			| Boolean						| true => Has Yearpicker						|
| **date**			| Boolean						| true => Has Datepicker						|
| **disabledDays**	| Array							| Array containing day indices (1-7) to disable	|
| **time**			| Boolean						| true => Has Timepicker						|
| **clearButton**	| Boolean						| true => Show Clear Button						|
| **nowButton**		| Boolean						| true => Show Now Button						|
| **switchOnClick**	| Boolean						| true => Switch view on click (default: false) |
| **cancelText**	| String						| Text for the cancel button (default: Cancel)	|
| **okText**		| String						| Text for the OK button (default: OK)			|
| **clearText**		| String						| Text for the Clear button (default: Clear)	|
| **nowText**		| String						| Text for the Now button (default: Now)		|
| **triggerEvent**		| String						| Event to fire the calendar (default: focus)		|
| **monthPicker**	| Boolean						| true => Act as monthpicker (hide calendar) (default: false) |



### Events

| Name				| Parameters				| Description										|
| ----------------- | ------------------------- | ------------------------------------------------- |
| **beforeChange**	| event, date				| OK button is clicked								|
| **change**		| event, date				| OK button is clicked and input value is changed	|
| **yearSelected**	        | event, date			        | New year is selected								|
| **dateSelected**	| event, date				| New date is selected								|
| **open**	        | event				        | datepicker is opened								|
| **close**	        | event				        | datepicker is closed								|


### Methods

        $('input').bootstrapMaterialDatePicker('setDate', moment());

| Name				| Parameter					| Description					|
| ----------------- | ------------------------- | ----------------------------- |
| **setDate**		| (String\|Date\|Moment)	| Set initial date				|
| **setMinDate**	| (String\|Date\|Moment)	| Set minimum selectable date	|
| **setMaxDate**	| (String\|Date\|Moment)	| Set maximum selectable date	|
| **destroy**		| NULL						| Destroy the datepicker		|
---
name: Bug report
about: Create a report to help us improve DAISY
title: ''
labels: ''
assignees: ''

---

**Describe the bug (current behaviour)**


**To reproduce**

1. 
2. 
3. 
4. 

**Expected behavior**


**Environment:**

 - Deployment: [e.g. daisy-test.lcsb.uni.lu]
 - Version [e.g. 1.1.0]
 - Browser [e.g. Firefox 68.x]
 
 
**Additional information (error logs, screenshots)**
---
name: Feature request
about: Suggest an new feature for DAISY
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**

**Describe the solution requested**


**Describe possible alternatives**


**Other info context**
