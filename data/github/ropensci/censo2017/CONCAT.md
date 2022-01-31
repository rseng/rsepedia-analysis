
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Censo 2017 (Paquete R)

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable-1)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![GH-actions](https://github.com/ropensci/censo2017/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/censo2017/actions)
[![codecov](https://codecov.io/gh/ropensci/censo2017/branch/main/graph/badge.svg?token=XI59cmGd15)](https://codecov.io/gh/ropensci/censo2017)
[![CRAN
status](https://www.r-pkg.org/badges/version/censo2017)](https://CRAN.R-project.org/package=censo2017)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4277761.svg)](https://doi.org/10.5281/zenodo.4277761)
[![Buy Me a
Coffee](https://img.shields.io/badge/buymeacoffee-pacha-yellow)](https://www.buymeacoffee.com/pacha?via=github)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/414_status.svg)](https://github.com/ropensci/software-review/issues/414)
<!-- badges: end -->

# Acerca de

Provee un acceso conveniente a mas de 17 millones de registros de la
base de datos del Censo 2017. Los datos fueron importados desde el DVD
oficial del INE usando el [Convertidor
REDATAM](https://github.com/discontinuos/redatam-converter/) creado por
Pablo De Grande y ademas se proporcionan los mapas que acompanian a
estos datos. Estos mismos datos en DVD posteriormente quedaron
disponibles en las [Bases de Datos del
INE](https://www.ine.cl/estadisticas/sociales/censos-de-poblacion-y-vivienda/poblacion-y-vivienda).

Despues de la primera llamada a `library(censo2017)` se le pedira al
usuario que descargue la base usando `censo_descargar_base()` y se puede
modificar la ruta de descarga con la variable de entorno
`CENSO_BBDD_DIR`. La variable de entorno se puede crear con
`usethis::edit_r_environ()`.

La documentacion esta disponible en
<https://docs.ropensci.org/censo2017/>.

# Publico objetivo

Estudiantes, academicos e investigadores que necesiten un acceso
conveniente a datos censales directamente en R o RStudio.

# Requerimientos de instalacion

Este paquete necesita 3.5 GB libres para la crear la base de datos
localmente.

# Instalacion

Version estable

    install.packages("censo2017")

Version de desarrollo

    # install.packages("remotes")
    remotes::install_github("ropensci/censo2017")

# Valor agregado sobre los archivos SHP y REDATAM del INE

Esta version de la base de datos del Censo 2017 presenta algunas
diferencias respecto de la original que se obtiene en DVD y corresponde
a una version DuckDB derivada a partir de los Microdatos del Censo 2017
en formato DVD.

La modificacion sobre los archivos originales, que incluyen geometrias
detalladas disponibles en [Cartografias
Censo2017](https://github.com/ropensci/censo2017-cartografias),
consistio en unir todos los archivos SHP regionales en una unica tabla
por nivel (e.g en lugar de proveer `R01_mapa_comunas`, …,
`R15_mapa_comunas` combine las 15 regiones en una unica tabla
`mapa_comunas`).

Los cambios concretos respecto de la base original son los siguientes:

-   Nombres de columna en formato “tidy” (e.g. `comuna_ref_id` en lugar
    de `COMUNA_REF_ID`).
-   Agregue los nombres de las unidades geograficas (e.g. se incluye
    `nom_comuna` en la tabla `comunas` para facilitar los filtros).
-   Aniadi la variable `geocodigo` a la tabla de `zonas`. Esto facilita
    mucho las uniones con las tablas de mapas en SQL.
-   Tambien inclui las observaciones 16054 to 16060 en la variable
    `zonaloc_ref_id`. Esto se debio a que era necesario para crear una
    llave foranea desde la tabla `mapa_zonas` (ver repositorio
    [Cartografias
    Censo2017](https://github.com/ropensci/censo2017-cartografias)) y
    vincular el `geocodigo` (no todas las zonas del mapa estan presentes
    en los datos del Censo).

Ademas de los datos del Censo, inclui la descripcion de las variables en
formato tabla (y no en XML como se obtiene del DVD). La ventaja de esto
es poder consultar rapidamente lo que significan los codigos de
variables y su etiquetado, por ejemplo:

``` r
# con la bbdd instalada
library(censo2017)
library(dplyr)

censo_tabla("variables") %>% 
  filter(variable == "p01")
#> # A tibble: 1 x 5
#>   tabla     variable descripcion      tipo    rango 
#>   <chr>     <chr>    <chr>            <chr>   <chr> 
#> 1 viviendas p01      Tipo de Vivienda integer 1 - 10

censo_tabla("variables_codificacion") %>% 
  filter(variable == "p01")
#> # A tibble: 12 x 4
#>    tabla     variable valor descripcion                                         
#>    <chr>     <chr>    <int> <chr>                                               
#>  1 viviendas p01          1 Casa                                                
#>  2 viviendas p01          2 Departamento en Edificio                            
#>  3 viviendas p01          3 Vivienda Tradicional Indígena (Ruka, Pae Pae u Otra…
#>  4 viviendas p01          4 Pieza en Casa Antigua o en Conventillo              
#>  5 viviendas p01          5 Mediagua, Mejora, Rancho o Choza                    
#>  6 viviendas p01          6 Móvil (Carpa, Casa Rodante o Similar)               
#>  7 viviendas p01          7 Otro Tipo de Vivienda Particular                    
#>  8 viviendas p01          8 Vivienda Colectiva                                  
#>  9 viviendas p01          9 Operativo Personas en Tránsito (No Es Vivienda)     
#> 10 viviendas p01         10 Operativo Calle (No Es Vivienda)                    
#> 11 viviendas p01         11 Valor Perdido                                       
#> 12 viviendas p01          0 No Aplica
```

# Relacion de Censo 2017 con Chilemapas

Todos los datos de estos repositorios contemplan 15 regiones pues los
archivos del Censo se entregan de esta forma y este paquete esta 100%
orientado a facilitar el acceso a datos.

Por su parte, [chilemapas](https://docs.ropensci.org/censo2017) se
centra unicamente en los mapas y tambien usa las cartografias del DVD
del Censo para entregar mapas simplificados (de menor detalle y mas
livianos). Chilemapas cuenta con una transformacion de codigos para dar
cuenta de la creacion de la Region de Niuble.

En resumen, censo2017 permite construir estadisticas demograficas y
chilemapas ayuda a mostrarlas en un mapa usando ggplot2 (u otro paquete
como tmap).

# Cita este trabajo

Si usas `censo2017` en trabajos academicos u otro tipo de publicacion
por favor usa la siguiente cita:

    Mauricio Vargas (2020). censo2017: Base de Datos de Facil Acceso del Censo
      2017 de Chile (2017 Chilean Census Easy Access Database). R package version
      0.1. https://docs.ropensci.org/censo2017/

Entrada para BibTeX:

    @Manual{,
      title = {censo2017: Base de Datos de F\'acil Acceso del Censo 2017 de Chile
    (2017 Chilean Census Easy Access Database)},
      author = {Mauricio Vargas},
      year = {2020},
      note = {R package version 0.1},
      url = {https://docs.ropensci.org/censo2017/},
      doi = {10.5281/zenodo.4277761}
    }

# Contribuciones

Para contribuir a este proyecto debes estar de acuerdo con el [Codigo de
Conducta de rOpenSci](https://ropensci.org/code-of-conduct/). Me es util
contar con mas ejemplos, mejoras a las funciones y todo lo que ayude a
la comunidad. Si tienes algo que aportar me puedes dejar un issue o pull
request.

# Agradecimientos

Muchas gracias a Juan Correa por su asesoria como geografo experto.

# Aportes

Si quieres donar para aportar al desarrollo de este y mas paquetes Open
Source, puedes hacerlo en [Buy Me a
Coffee](https://www.buymeacoffee.com/pacha/).
# Version 0.5

- Adds all of rOpenSci's suggestion, which means great improvements such as
 better documentation, consistent syntax, does not depend on DuckDB version, etc.
 See https://github.com/ropensci/software-review/issues/414 for the full detail.
 
# Version 0.4

- Works with duckdb 0.3.4
- Removes the databases used with older censo2017 versions

# Version 0.3

- Moves local database location according to CRAN request
- Requires R 4.0
- Uses DuckDB instead of SQLite

# Version 0.2

- Adds `vignettes/` to `.Rbuildignore`.
- Complies with CRAN policies regarding Suggests.
- Vignettes are now available from gh-pages.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# Censo 2017 (Paquete R)

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable-1)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![GH-actions](https://github.com/ropensci/censo2017/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/censo2017/actions)
[![codecov](https://codecov.io/gh/ropensci/censo2017/branch/main/graph/badge.svg?token=XI59cmGd15)](https://codecov.io/gh/ropensci/censo2017)
[![CRAN status](https://www.r-pkg.org/badges/version/censo2017)](https://CRAN.R-project.org/package=censo2017)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4277761.svg)](https://doi.org/10.5281/zenodo.4277761)
[![Buy Me a Coffee](https://img.shields.io/badge/buymeacoffee-pacha-yellow)](https://www.buymeacoffee.com/pacha?via=github)
[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/414_status.svg)](https://github.com/ropensci/software-review/issues/414)
<!-- badges: end -->

# Acerca de

Provee un acceso conveniente a mas de 17 millones de registros de la base de datos del Censo 2017. Los datos fueron importados desde el DVD oficial del INE usando el [Convertidor REDATAM](https://github.com/discontinuos/redatam-converter/) creado por Pablo De Grande y ademas se proporcionan los mapas que acompanian a estos datos. Estos mismos datos en DVD posteriormente quedaron disponibles en las [Bases de Datos del INE](https://www.ine.cl/estadisticas/sociales/censos-de-poblacion-y-vivienda/poblacion-y-vivienda).

Despues de la primera llamada a `library(censo2017)` se le pedira al usuario que descargue la base usando `censo_descargar_base()` y se puede modificar la ruta de descarga con la variable de entorno `CENSO2017_DIR`. La variable de entorno se puede crear con `usethis::edit_r_environ()`.

La documentacion esta disponible en https://docs.ropensci.org/censo2017/.

# Publico objetivo

Estudiantes, academicos e investigadores que necesiten un acceso conveniente a datos censales directamente en R o RStudio.

# Requerimientos de instalacion

Este paquete necesita 3.5 GB libres para la crear la base de datos localmente.

# Instalacion

Version estable
```
install.packages("censo2017")
```

Version de desarrollo
```
# install.packages("remotes")
remotes::install_github("ropensci/censo2017")
```

# Valor agregado sobre los archivos SHP y REDATAM del INE

Esta version de la base de datos del Censo 2017 presenta algunas diferencias respecto de la original que se obtiene en DVD y corresponde a una version DuckDB derivada a partir de los Microdatos del Censo 2017 en formato DVD. 

La modificacion sobre los archivos originales, que incluyen geometrias detalladas disponibles en [Cartografias Censo2017](https://github.com/ropensci/censo2017-cartografias), consistio en unir todos los archivos SHP regionales en una unica tabla por nivel (e.g en lugar de proveer `R01_mapa_comunas`, ..., `R15_mapa_comunas` combine las 15 regiones en una unica tabla `mapa_comunas`).

Los cambios concretos respecto de la base original son los siguientes:

* Nombres de columna en formato "tidy" (e.g. `comuna_ref_id` en lugar de `COMUNA_REF_ID`).
* Agregue los nombres de las unidades geograficas (e.g. se incluye `nom_comuna` en la tabla `comunas` para facilitar los filtros).
* Aniadi la variable `geocodigo` a la tabla de `zonas`. Esto facilita mucho las uniones con las tablas de mapas en SQL.
* Tambien inclui las observaciones 16054 to 16060 en la variable `zonaloc_ref_id`. Esto se debio a que era necesario para crear una llave foranea desde la tabla `mapa_zonas` (ver repositorio [Cartografias Censo2017](https://github.com/ropensci/censo2017-cartografias)) y vincular el `geocodigo` (no todas las zonas del mapa estan presentes en los datos del Censo).

Ademas de los datos del Censo, inclui la descripcion de las variables en formato tabla (y no en XML como se obtiene del DVD). La ventaja de esto es poder consultar rapidamente lo que significan los codigos de variables y su etiquetado, por ejemplo:
```{r message=FALSE, warning=FALSE}
# con la bbdd instalada
library(censo2017)
library(dplyr)

censo_tabla("variables") %>% 
  filter(variable == "p01")

censo_tabla("variables_codificacion") %>% 
  filter(variable == "p01")
```

# Relacion de Censo 2017 con Chilemapas

Todos los datos de estos repositorios contemplan 15 regiones pues los archivos del Censo se entregan de esta forma y este paquete esta 100% orientado a facilitar el acceso a datos.

Por su parte, [chilemapas](https://docs.ropensci.org/censo2017) se centra unicamente en los mapas y tambien usa las cartografias del DVD del Censo para entregar mapas simplificados (de menor detalle y mas livianos). Chilemapas cuenta con una transformacion de codigos para dar cuenta de la creacion de la Region de Niuble.

En resumen, censo2017 permite construir estadisticas demograficas y chilemapas ayuda a mostrarlas en un mapa usando ggplot2 (u otro paquete como tmap).

# Cita este trabajo

Si usas `censo2017` en trabajos academicos u otro tipo de publicacion por favor usa la siguiente cita:

```
Mauricio Vargas (2020). censo2017: Base de Datos de Facil Acceso del Censo
  2017 de Chile (2017 Chilean Census Easy Access Database). R package version
  0.1. https://docs.ropensci.org/censo2017/
```

Entrada para BibTeX:

```
@Manual{,
  title = {censo2017: Base de Datos de F\'acil Acceso del Censo 2017 de Chile
(2017 Chilean Census Easy Access Database)},
  author = {Mauricio Vargas},
  year = {2020},
  note = {R package version 0.1},
  url = {https://docs.ropensci.org/censo2017/},
  doi = {10.5281/zenodo.4277761}
}
```

# Contribuciones

Para contribuir a este proyecto debes estar de acuerdo con el [Codigo de Conducta de rOpenSci](https://ropensci.org/code-of-conduct/). Me es util contar con mas ejemplos, mejoras a las funciones y todo lo que ayude a la comunidad. Si tienes algo que aportar me puedes dejar un issue o pull request.

# Agradecimientos

Muchas gracias a Juan Correa por su asesoria como geografo experto.

# Aportes

Si quieres donar para aportar al desarrollo de este y mas paquetes Open Source, puedes hacerlo en [Buy Me a Coffee](https://www.buymeacoffee.com/pacha/).
---
title: "Uso basico del paquete censo2017"
author: "Mauricio Vargas S."
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Uso basico del paquete censo2017}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  cache = FALSE,
  collapse = TRUE,
  eval = TRUE,
  comment = "#>"
)
```

# Introduccion

Este paquete se integra perfectamente con el tidyverse, con el que daremos un ejemplo muy basico para mostrar sus principales funciones.

# Aproximacion de la poblacion con el grado de doctor en la Region del Bio Bio.

Se procedera a obtener una aproximacion usando dplyr ya que, puede haber personas que no sean de la comuna y aparezcan censadas en ella. Sin embargo, cabe aclarar que para esta ocasion no haremos el filtro que corrige esto. Nuestra idea es mantener el ejemplo lo mas simple posible.

Primero que todo, se cargan los paquetes necesarios.

• censo2017: Proporciona los datos censales para poder generar las tablas y graficos de este ejemplo,
• dplyr: Facilta filtrar datos en una tabla, unir distintas tablas y en general todas las tareas de limpieza y transformacion de datos. 
* ggplot2: Nos permite graficar usando el concepto de "gramatica de graficos", es decir que podemos ir creando graficos incrementales y controlar los ejes, titulos y demas elementos por separado.
• chilemapas: Nos entrega mapas terrestres con topologias simplificadas. Esto significa, que este paquete al contener poligonos que generan graficos, y no tablas, nos permite un uso mas sencillo para lo que queremos realizar, que es darle un complemento visual a la información cargada desde censo2017.

```{r, warning=FALSE, message=FALSE}
library(censo2017)
library(dplyr)
library(ggplot2)
library(chilemapas)
```

Hay que realizar algunos cruces de tablas, de manera de filtrar la region que nos interesa.

Comenzamos con la tabla zonas: generamos la provincia a partir del geocodigo y luego filtro para unir hasta llegar a la tabla personas. Nos interesa utilizar la variable `p15`, cuya descripcion esta en la tabla `variables` y cuya codificacion aparece en la tabla `variables_codificacion`.

```{r, warning=FALSE, message=FALSE, eval=FALSE}
nivel_educacional_biobio <- tbl(censo_conectar(), "zonas") %>% 
  mutate(
    region = substr(as.character(geocodigo), 1, 2),
    comuna = substr(as.character(geocodigo), 1, 5)
  ) %>% 
  filter(region == "08") %>% 
  select(comuna, geocodigo, zonaloc_ref_id) %>%
  inner_join(select(tbl(censo_conectar(), "viviendas"), zonaloc_ref_id, vivienda_ref_id), by = "zonaloc_ref_id") %>%
  inner_join(select(tbl(censo_conectar(), "hogares"), vivienda_ref_id, hogar_ref_id), by = "vivienda_ref_id") %>%
  inner_join(select(tbl(censo_conectar(), "personas"), hogar_ref_id, nivel_educ = p15), by = "hogar_ref_id") %>%
  collect()
```

Con lo anterior, los niveles educacionales de las personas censadas se pueden agrupar por comuna y obtener la cuenta proporcionada en base a la suma total.
```{r, warning=FALSE, message=FALSE, eval=FALSE}
nivel_educacional_biobio <- nivel_educacional_biobio %>% 
  group_by(comuna, nivel_educ) %>%
  summarise(cuenta = n()) %>%
  group_by(comuna) %>%
  mutate(proporcion = cuenta / sum(cuenta))
```

Vemos los datos antes de continuar.
```{r}
nivel_educacional_biobio
```

Creamos la variable mapa_biobio haciendo un filtro para obtener unicamente los datos de la region con codigo "08" (region del Bio Bio) desde mapa_comunas. Luego de eso, haremos un left join (union de tablas manteniendo todas las filas de la tabla izquierda o inicial) desde la tabla chilemapas, donde obtendremos el mapa de la provincia, y la uniremos con los datos coincidentes segon el codigo_comuna de la tabla censo2017.
```{r, warning=FALSE, message=FALSE}
mapa_biobio <- mapa_comunas %>% 
  filter(codigo_region == "08") %>% 
  left_join(nivel_educacional_biobio, by = c("codigo_comuna" = "comuna"))
```

Ahora que cargamos toda la informacion necesaria en R desde la base de datos, debemos cerrar la conexion SQL (importante).
```{r, warning=FALSE, message=FALSE}
censo_desconectar()
```

Finalmente procedemos a generar el mapa.

Primero, creamos la variable colors, en el que incluiremos los codigos hexadecimales de los colores que utilizaremos al momento de crear el mapa.

Luego de hecho esto, utilizamos geom_sf del paquete ggplot2, que se usa para visualizar objetos de caracteristicas simples (sf = simple features). Geom_sf dibujara diferentes objetos geometricos a partir de una columna de tipo 'sf' que debe estar presente en los datos.

Seleccionamos el codigo_comuna y geometry (que contiene los poligonos que componen cada region) desde el mapa_biobio, que creamos anteriormente. Volvemos a hacer un left_join ahora de mapa_biobio seleccionando unicamente codigo_comuna, nivel_educ y proporcion.
```{r, fig.width=10, warning=FALSE, message=FALSE}
colors <- c("#DCA761","#C6C16D","#8B9C94","#628CA5","#5A6C7A")

g <- ggplot() +
  geom_sf(data = mapa_biobio %>% 
            select(codigo_comuna, geometry) %>% 
            left_join(
              mapa_biobio %>% 
                filter(nivel_educ == 14) %>% 
                select(codigo_comuna, nivel_educ, proporcion),
              by = "codigo_comuna"
            ),
          aes(fill = proporcion, geometry = geometry),
          size = 0.1) +
  scale_fill_gradientn(colours = rev(colors), name = "Porcentaje") +
  labs(title = "Porcentaje de habitantes con el grado de doctor\npor comuna en la Region del Bio Bio") +
  theme_minimal(base_size = 13)

g
```

Notas:

* El uso de `tbl()` y `collect()` en la primera parte se podra entender mejor leyendo, por ejemplo, [A Crash Course on PostgreSQL for R Users](https://pacha.dev/blog/2020/08/09/a-crash-course-on-postgresql-for-r-users/).
* En la segunda parte se usa `censo_tabla()` ya que SQL almacena la columna `geometry` (de tipo poligono) como cadena de texto mientras que R lee poligonos sin problema.
* En la tercera parte hago un join entre el mapa completo y la tabla con quienes tienen el grado de doctor. Este paso, aunque pueda parecer redundante, es necesario si quiero mostrar las zonas con 0 doctores y si lo omito se borran algunas zonas del mapa.
* El mapa que se genero usando las funciones de chilemapas podria haber generado con las cartografias oficiales del Censo (ver repositorio cartografias-censo2017. Esta alternativa entrega un mayor nivel de detalle, pero requiere mayor esfuerzo para leer las cartografias y el tiempo requerido para generar los mapas aumenta fuertemente.

# Ejercicios para el usuario

1. Realizar un grafico similar al del ejemplo pero a nivel de zona censal.
2. Explorar la columna `p10` en la tabla `personas` y realizar un grafico que de cuenta de la poblacion efectiva de la comuna (e.g. mejorando el problema de personas que podrian no ser de la comuna en el ejemplo).
3. Agregar datos al mapa sin usar `chilemapas`. Una forma de hacerlo es la siguiente:

```{r, warning=FALSE, message=FALSE, eval=FALSE}
mapa_biobio <- censo_tabla("mapa_comunas") %>%
  filter(region == "08") %>% 
  left_join(nivel_educacional_biobio, by = "comuna")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/censo2017-package.R
\docType{package}
\name{censo2017-package}
\alias{censo2017}
\alias{censo2017-package}
\title{censo2017: Base de Datos de Facil Acceso del Censo 2017 de Chile
    (2017 Chilean Census Easy Access Database)}
\description{
Provee un acceso conveniente a mas de 17 millones de registros
    de la base de datos del Censo 2017. Los datos fueron importados desde
    el DVD oficial del INE usando el Convertidor REDATAM creado por Pablo De
    Grande. Esta paquete esta documentado intencionalmente en castellano
    asciificado para que funcione sin problema en diferentes plataformas.
    (Provides convenient access to more than 17 million records from the
    Chilean Census 2017 database. The datasets were imported from the official
    DVD provided by the Chilean National Bureau of Statistics by using the
    REDATAM converter created by Pablo De Grande and in addition it includes the
    maps accompanying these datasets.)
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/censo2017/}
  \item Report bugs at \url{https://github.com/ropensci/censo2017/issues/}
}

}
\author{
\strong{Maintainer}: Mauricio Vargas \email{mavargas11@uc.cl} (\href{https://orcid.org/0000-0003-1017-7574}{ORCID})

Other contributors:
\itemize{
  \item Juan Correa [contributor]
  \item Maria Paula Caldas (rOpenSci) [reviewer]
  \item Frans van Dunee (rOpenSci) [reviewer]
  \item Melina Vidoni (rOpenSci) [reviewer]
  \item Constanza Manriquez (revision independiente de las vinietas) [reviewer]
  \item  Instituto Nacional de Estadisticas (INE) [data contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connect.R
\name{censo_conectar}
\alias{censo_conectar}
\title{Conexion a la Base de Datos del Censo}
\usage{
censo_conectar(dir = censo_path())
}
\arguments{
\item{dir}{La ubicacion de la base de datos en el disco. Por defecto es
\code{censo2017} en la carpeta de datos del usuario de R o la variable de entorno
\code{CENSO2017_DIR} si el usuario la especifica.}
}
\description{
Devuelve una conexion a la base de datos local. Esto corresponde a una
conexion a una base DuckDB compatible con DBI. A diferencia de
\code{\link[=censo_tabla]{censo_tabla()}}, esta funcion es mas flexible y se puede usar con
dbplyr para leer unicamente lo que se necesita o directamente con DBI para
usar comandos SQL.
}
\examples{
\dontrun{
 DBI::dbListTables(censo_conectar())

 DBI::dbGetQuery(
  censo_conectar(),
  'SELECT * FROM comunas WHERE provincia_ref_id = 1'
 )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connect.R
\name{censo_tabla}
\alias{censo_tabla}
\title{Tablas Completas de la Base de Datos del Censo}
\usage{
censo_tabla(tabla)
}
\arguments{
\item{tabla}{Una cadena de texto indicando la tabla a extraer}
}
\value{
Un tibble
}
\description{
Devuelve una tabla completa de la base de datos. Para entregar datos
filtrados previamente se debe usar \code{\link[=censo_conectar]{censo_conectar()}}.
}
\examples{
\dontrun{ censo_tabla("comunas") }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{censo_descargar}
\alias{censo_descargar}
\title{Descarga la Base de Datos del Censo a tu Computador}
\usage{
censo_descargar(ver = NULL)
}
\arguments{
\item{ver}{La version a descargar. Por defecto es la ultima version
disponible en GitHub. Se pueden ver todas las versiones en
\url{https://github.com/pachamaltese/censo2017/releases}.}
}
\description{
Este comando descarga la base de datos completa como un unico archivo zip que
se descomprime para crear la base de datos local. Si no quieres descargar la
base de datos en tu home, ejecuta usethis::edit_r_environ() para crear la
variable de entorno CENSO2017_DIR con la ruta.
}
\examples{
\dontrun{ censo_descargar() }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connect.R
\name{censo_desconectar}
\alias{censo_desconectar}
\title{Desconecta la Base de Datos del Censo}
\usage{
censo_desconectar()
}
\description{
Una funcion auxiliar para desconectarse de la base de datos.
}
\examples{
censo_desconectar()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/censo2017-package.R
\docType{data}
\name{nivel_educacional_biobio}
\alias{nivel_educacional_biobio}
\title{Poblacion por Nivel Educacional en la Region del Bio Bio}
\format{
Un tibble con 860 observaciones en las siguientes 4 variables.
\describe{
\item{\code{comuna}}{codigo de comuna (15 regiones)}
\item{\code{nivel_educ}}{maximo nivel educacional alcanzado (ver la vinieta
con los links a la descripcion de codigos)}
\item{\code{cuenta}}{cantidad de personas censadas en la comuna}
\item{\code{proporcion}}{porcentaje que representan las personas censadas en
la comuna}
}
}
\description{
Proporciona la cuenta y porcentaje por comuna de las personas de
la Region del Bio Bio de acuerdo al maximo nivel educacional que reportan
(e.g. primaria, secundaria, universitaria, etc.)
}
\examples{
nivel_educacional_biobio

\dontrun{
# replicar el resultado usando dplyr directamente con SQL
# es ligeramente distinto a las vinietas que explican esta misma tabla
nivel_educacional_biobio <- tbl(censo_conectar(), "zonas") \%>\% 
 mutate(
  region = substr(as.character(geocodigo), 1, 2),
  comuna = substr(as.character(geocodigo), 1, 5)
 ) \%>\% 
 filter(region == "08") \%>\% 
 select(comuna, geocodigo, zonaloc_ref_id) \%>\%
 inner_join(select(tbl(censo_conectar(), "viviendas"),
  zonaloc_ref_id, vivienda_ref_id), by = "zonaloc_ref_id") \%>\%
 inner_join(select(tbl(censo_conectar(), "hogares"),
  vivienda_ref_id, hogar_ref_id), by = "vivienda_ref_id") \%>\%
 inner_join(select(tbl(censo_conectar(), "personas"),
  hogar_ref_id, nivel_educ = p15), by = "hogar_ref_id") \%>\%
 group_by(comuna, nivel_educ) \%>\%
 summarise(cuenta = n()) \%>\%
 group_by(comuna) \%>\%
 mutate(proporcion = cuenta * (1 / sum(cuenta))) \%>\% 
 arrange(comuna, nivel_educ)}
}
\author{
Elaboracion propia con base en datos desagregados del Censo
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove.R
\name{censo_eliminar}
\alias{censo_eliminar}
\title{Elimina la Base de Datos del Censo de tu Computador}
\usage{
censo_eliminar(preguntar = TRUE)
}
\arguments{
\item{preguntar}{Si acaso se despliega un menu para confirmar la accion de
borrar cualquier base del censo existente. Por defecto es verdadero.}
}
\description{
Elimina el directorio \code{censo2017} y todos sus contenidos, incluyendo versiones
de la base de datos del Censo creadas con cualquier version de 'DuckDB'.
}
\examples{
\dontrun{ censo_eliminar() }
}
