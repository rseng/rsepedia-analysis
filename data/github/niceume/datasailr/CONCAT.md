# About DataSailr

DataSailr is an R package which enables row by row data manipulation. 

The data manipulation process is instructed using DataSailr script, which is designed specifically for data manipulation. The internal data manipulation engine is a library called libsailr, which is written in C/C++. If you are interested in DataSailr's internal implementation, please also have a look at libsailr library.

This package is still in development. Please let me know if you find bugs or inconsistent behaviors. Please also read community guidelines at the end of this document, when you wish to contribute to this software, report issues or problems, or need supports.


## Motivation

Data processing and data cleaning are essential steps before applying statistical or machine learning procedures. R provides a flexible way for data processing using vectors. Additional R packages also provide other ways for manipulating data such as using SQL and using chained functions. The datasailr package provides yet another way to process data in a row by row manner using data manipulation oriented script, DataSailr script. The script enables users to instruct data processing in a row wise manner, and the script is intuitive and easy to learn. This row by row approach is different from R's built-in functionality, which provides data processing through column vector calculations. Row wise approach does not require vectors, and can be intuitive compared to R's default column vector approach.

DataSailr script is a domain specific language for data processing. Using domain specific language for data manipulation is a similar approach to SAS software (SAS and all other SAS Institute Inc. product or service names are registered trademarks or trademarks of SAS Institute Inc. in the USA and other countries). SAS software provides DATA blocks, within which users can write scripts that are specific for data manipulation in a row wise manner. The separation of data manipulation steps and statistical procedures has a benefit of improved readability. Its row wise data manipulation script is easy to understand and learn. DataSailr brings the same kind of experience to R.

As mentioned, there are R packages that have been improving data manipulation, such as sqldf and dplyr. The sqldf package enables users to write SQL for data manipulation. The dplyr package enables users to write data manipulation procedures in a sequential way by chaining functions without thinking much about column vectors. DataSailr enables the same kind of things with a single DataSailr script.


## How to install

### CRAN

DataSailr package is available on CRAN.

```
# Start R interpreter
install.packages("datasailr")
```


### Binaries

Binaries are provided at "https://datasailr.io/download"


### Install from source

For compilation, autotools, bison, flex and some other UNIX programs are required.

DataSailr_toolchain repository provides a way to create CRAN package, and to install it. CRAN package requires some rules. For example, std::cout is not used, but Rcpp::Rcout needs to be used. Also, DataSailr repository itself does not include libsailr and Onigmo. When including those libraries for CRAN package, additional copy right holders need to be listed in DESCRIPTION file. If you do not need these modifications, installing from DataSailr repository itself is also an acceptable way..


* Install from datasailr_toolchain repository

```
git clone https://github.com/niceume/datasailr_toolchain.git
cd datasailr_toolchain
./create_datasailr_pkg.sh
./create_packages.sh
```

* Install from datasailr repository

```
git clone <datasailr repository>
cd datasailr/src
./update_3rd_party_library.sh  # This downloads Onigmo library.
cd ..
autoconf  # Generate configure script.
cd ..

R CMD INSTALL datasailr --preclean --no-multiarch --build
```


## How to run DataSailr script

datasailr::sail() is the main function of this package. This function takes data.frame as the 1st argument and DataSailr script as the 2nd argument. (Note that the 1st argument is data.frame, so you can use this package with %>% operator (in magrittr package), and combine this with functions of dplyr.)

```
library(datasailr)

datasailr::sail( df, '
<DataSailr script>
'
)
```


## Quick introduction

In this quick introduction, built-in R dataframe called mtcars is going to be used. The mtcars datafram has automobile information in 1970's. The detail of this dataframe is described here. (https://www.rdocumentation.org/packages/datasets/versions/3.6.1/topics/mtcars)

* Built-in dataframe, mtcars.

```
data(mtcars)
head(mtcars,2)
```

(OUTPUT)

|              | mpg | cyl | disp | hp  | drat |   wt  | qsec  |vs  | am |gear |carb
---------------|-----|-----|------|-----|------|-------|-------|----|----|-----|-----
|Mazda RX4     | 21  | 6   | 160  | 110 | 3.9  | 2.620 | 16.46 | 0  | 1  | 4   | 4
|Mazda RX4 Wag | 21  | 6   | 160  | 110 | 3.9  | 2.875 | 17.02 | 0  | 1  | 4   | 4


### Code example

Before executing the following code, install DataSailr package (and magrittr package.) Installation instruction is mentioned later in this README (Binary package is the easiest way to try, until this package becomes available on CRAN).

* (Example)

```
library(datasailr)
library(magrittr)

mtcars %>% datasailr::sail('
if( hp > 145 ){  power = "high" }
else if( 145 >= hp && hp > 0){ power = "low" }
else{  print("hp variable has missing value")}

germany = re/(^Merc|^Porsche|^Volvo)/
usa = re/(^Hornet|^Cadillac|^Lincoln|^Chrysler|^Dodge|^AMC|^Camaro|^Chevrolet|^Pontiac|^Ford)/
japan = re/(^Mazda|^Datsun|^Honda|^Toyota)/

if ( _rowname_ =~ germany ) { country = "Germany" ; type = rexp_matched(1); }
else if( _rowname_ =~ usa ) { country = "USA"  ; type = rexp_matched(1);  }
else if( _rowname_ =~ japan ) { country = "Japan"  ; type = rexp_matched(1); }
else { country = "Other" }
') %>% head(10)
```


* (OUTPUT)

	
|   | power | country | type
----|-------|---------|-----
| 1 | low   | Japan   | Mazda	
| 2 | low   | Japan   | Mazda	
| 3 | low   | Japan   | Datsun	
| 4 | low   | USA     | Hornet	
| 5 | high  | USA     | Hornet	
| 6 | low   | Other   |
| 7 | high  | Other   |
| 8 | low   | Germany | Merc	
| 9 | low   | Germany | Merc	
| 10 | low  | Germany | Merc


### Example code explanation

* library(datasailr)
    + Make DataSailr package available.
* library(magrittr)
    + This is not required by DataSailr. In the above code, %>% , pipe operator, is used, and it's defined in magrittr package.
* mtcars %>% datasailr::sail(' .... 
    + The 1st argument of datasailr::sail is dataframe. In this case, %>% passes dataframe to the 1st argument.
    + The 2nd argument of datasailr::sail is DataSailr code.
* if( hp > 145 ){ ....
    + if-else statement. 
* power = "high" 
    + assign "high" to power column.
* germany = re/(^Merc|^Porsche|^Volvo)/
    + Regular expression literal is re/pattern/ .
* if ( \_rowname\_ =~ germany ){ ....
    + string =~ regular expression executes regular expression matching.
    + \_rowname\_ column is one of automatic variables. This holds the row names of the original R dataframe.
* rexp_matched(1)
    + Calling function of rexp_matched()
    + Extract the 1st matching group from the last regular expression matching.



## Grammar of Sailr/DataSailr script

Before describing grammar details, I need to mention that DataSailr script and Sailr script are essentially the same script. DataSailr script, script for DataSailr, is just an extended version of Sailr script(, script of libsailr). DataSailr and Sailr have the same grammar, and DataSailr just has additional built-in functions, such as push!() and discard!(). What is mentioned about Sailr script can be completely applied to DataSailr script.


Sailr script is especially intended for use of data manipulation. From the programming language view, the functionality is not enough and is not a general-purpose programming language, but you can write instructions and manipulate data in an intuitive way. For SAS(*) software users, the grammar of Sailr script is relatively easy to understand. Differences are mentioned later in this README.


(*) SAS is the world's largest privately held software company. SAS provides SAS software which is used by professional data analysts.


### Variables & Assignment operation [IMPORTANT]

* Variables in Sailr script and variables in statistics

In general-purpose programming language, variables usually mean identifiers to point to some memory or some objects. However in statistics, dataset column names are usually called variables. Sailr script follows the statistics way because it is designed not for general-purpose programming language, but is designed for processing dataset. Variables in codes correspond to the column with the same name. 

(e.g.) The following code print the value of column X.

```
print(X)
```

(e.g.) The following code assign the value of column X to column Y.

```
Y = X
```


### Types

Sailr script can deal with the following types.

1. Int
2. Double
    + Int and Double are interchangeable when needed.
3. String
    + String literal 
4. Boolean 
    + Only available within if-else conditions
5. Regular expression
    + Currently regular expressions are accepted only in literal. There is no such way to dynamically generate regular expression object.


R can handle vectors, list, S4 class objects, and so on. Among those R objects, dataframe usually consists of only vectors. The following list shows which data types are supported by Sailr. Supporting vectors of numbers and strings is enough in most of the cases of data analysis. 

* logical vector : Sailr imports TRUE as 1, and FALSE as 0.
* integer vector : Sailr deals with it as integer.
* numeric vector : Sailr deals with it as double.
* complex vector : Not supported
* character vector : Sailr deals with it as string.
* factor : Sailr deals with it as string.
    + Precisely, factor is a type called S3 class. S3 class is an R type with "class" attribute and some additional attributes (such as levels) required by the class. Internally factor is just an integer vector with these attributes. Attribute of "levels" is assumed to be CharacterVector in the current implementation.
* list : not supported.
* S3 class objects : As mentioned in factor, the base type of S3 is a basic R type. If the basic R type is supported by Sailr, the S3 type can be dealt from Sailr. However, Sailr's interpretation is not as the user's intention, but follows the internal data type. (i.e.) Factor is recognized as an integer vector.
* S4 class objects : Not supported
    + S4 objects are not supported. Those objects need to be converted to Sailr types beforehand. Moreover, S4 objects are only stored in list, and list will not be supported in any case.


| Column of R data.frame |  Sailr    | 
-------------------------|-----------------
| Integer Vector         | Int           |
| Double (Real) Vector   | Double        | 
| Character Vector       | String        |
| Boolean Vector         | Int (0/1)     |
| Factor                 | String        |
| List (No plan)         | No plan to be supported |


### Arithmetic operators


* Supported

The following operators can be applied to numbers (int + double).

Only plus (+) operator can deal with strings. It can concatenate strings. 

```
+    add
-    subtract
*    multiply
/    divide
**   power
^    power
=~   regular expression matching (  )
```

* Not supported

The following operators are not supported in Sailr

```
++  // Increment operator
--  // Decrement operator
+=
-=
*=
/=
```


### if-else statement

The following syntax is accepted. The former one is more familiar to programmers, though the latter is easier to read and may be better for data scientists.

```
if(bmi >= 40){
    weight_level = .
}else if( bmi >= 35 ){ 
    weight_level = .
}else if( bmi >= 30 ){
    weight_level = 3 
}else if( bmi >= 25 ){
    weight_level = 2
}else if( bmi >= 20 ){
    weight_level = 1
}else{
    weight_level = .
}
```

```
if(bmi >= 40){ weight_level = . } 
else if( bmi >= 35 ){ weight_level = . } 
else if( bmi >= 30 ){ weight_level = 3 } 
else if( bmi >= 25 ){ weight_level = 2 }
else if( bmi >= 20 ){ weight_level = 1 }
else { weight_level = . }
```


### Built-in functions

The following is the list of built-in functions.

```
print( str1, str2 ... )
num_to_str( num ) 
str_strip( str )
str_lstrip( str ) 
str_rstrip( str )
str_concat( str1, str2 ... ) // numbers are also accepted as input.
str_repeat( str , num )
str_subset( str, num, num )  // index starts from one. UTF8 compatible.
str_to_num( str )
rexp_matched( num )  // back reference
date_ymd( year, month, day )  // convert a combination of year, month and day into numbers that is a day from Jan 1st 1970.
date_ym_weekday_nth( year, month, weekday, nth )  // weekday should be "Sun", "Mon" ... 
date_add_n_years( unix_date, years )
date_add_n_months( unix_date, months )
date_add_n_days( unix_date, days )
date_format( unix_date, str_format )  // This format should follow C++'s std::chrono::format. "%m/%d/%y" or "%Y-%m-%d" is popular.
```

Currently, users cannot implement functions of Sailr script by themselves. 

If you want to implement your own function, you need to modify libsailr or libsailr's external function mechanism.

The libsailr C API is not mature. The interface for users is not organized. Please feel free to contact the author if you want your own function.



### Support for special numbers

#### Missing value support

For numbers (integer or double), single dot (.) represents missing value in code. Internally missing value is defined as -nan in C (which is generated by sqrt(-1)). See libsailr's node.c file and new_node_nan_double () function. ) Any calculation with this value always results in -nan.

For string, empty string ("") represents missing string.



#### Missing value (-nan) for if-condition

Note that in vanilla R, calculation with missing value returns missing value, and missing values act like TRUE in logical vector. In vanilla R, filtering of vectors and dataframes is done using logical vectors. If missing values are treated like TRUE, unintended elements or rows that have missing values are extracted. For example, you are calculating BMI (body mass index) from body weight and height. One person has a missing for body weight. His/her BMI will be calculated as missing value. You would like to categorize this person as "unknown BMI" category, then missingness should be checked first and these people should be categorized first. Otherwise, missing values in logical vector works like TRUE, so these people are categorized into some BMI categories, such as "high BMI" or "low BMI".

In Sailr, on the other hand, comparison operators with missing values return false except that missing is compared with missing.

```
library(datasailr)

col1 = c(1.7, 1.5 , 1.5 , 1.65 ,  1.5)
col2 = c(50 , 60  , 80  , NA   ,  100)
df = data.frame(height = col1, bw = col2)

code ='
bmi = bw / (height ^ 2)

if(bmi >= 40){ overweight = 5 } 
else if( bmi >= 35 ){ overweight = 4 } 
else if( bmi >= 30 ){ overweight = 3 } 
else if( bmi >= 25 ){ overweight = 2 }
else if( bmi >= 20 ){ overweight = 1 }
else if( bmi >= 0  ){ overweight = 0 }
else { overweight = . }
'

new_df = datasailr::sail( df, code , fullData = T);
new_df
```

* (OUTPUT)

```
  height  bw      bmi overweight
1   1.70  50 17.30104          0
2   1.50  60 26.66667          2
3   1.50  80 35.55556          4
4   1.65  NA       NA        NaN
5   1.50 100 44.44444          5
```


#### Infinite numbers

Infinite numbers are also supported. Internally this is implemented as Inf in C. The behavior is same as C.

(ref.) NAN and infinity in C Standard ( https://www.gnu.org/software/libc/manual/html_node/Infinity-and-NaN.html )
(ref.) This support seems to require at least C99. ( https://stackoverflow.com/questions/6235847/how-to-generate-nan-infinity-and-infinity-in-ansi-c )

```
x = 10
y = 0
z = -10

a = x / y /* INF */
b = z / y /* -INF */
c = (-1) ** (-1/2) /* sqrt(-1) => -nan */

print("a is ", a , "\n")
print("b is ", b , "\n")
print("c is ", c , "\n")

if( a >= 0 ){
  print("a is zero or more than zero. \n")
}else if( a < 0){
  print("a is less than zero. \n")
}else{
  print("a is not a number. \n")
}

if( b >= 0 ){
  print("b is zero or more than zero. \n")
}else if( b < 0){
  print("b is less than zero. \n")
}else{
  print("b is not a number. \n")
}

if( c >= 0 ){
  print("c is zero or more than zero. \n")
}else if(c < 0){
  print("c is less than zero. \n")
}else if(c == .){
  print("c is not a number. \n")
}
```

* OUTPUT

```
a is inf
b is -inf
c is -nan
a is zero or more than zero. 
b is less than zero. 
c is not a number. 
```


### UTF8 support and encodings

UTF8 encoding is strongly encouraged. The recommended way is to have dataframe as UTF8, and also have your code as UTF8. ASCII or latin1 is compatible with UTF8, so they must also work properly. (Your terminal encoding also need to be UTF8 for properly printing out UTF8 strings.)
 
In Sailr (libsailr), strings are stored just as a sequence of bytes. So even if you use other encodings, just copying and concatenating strings must work fine. Problems must happen when your source code and dataframe do not have the same encoding. Encoding of source code affects the string objects that are generated during parsing phase (=string objects that are generated before running code), dataframe encoding affects the string that are accessed via variables. 

For functions that process strings, UFF8-CPP library is used internally. ( UTF8-CPP:UTF-8 with C++ in a Portable Way https://github.com/nemtrif/utfcpp )

```
greeting = "Hello World"  // Only ascii string
str_len(greeting)  // 11
str_subset(1,5)  // Hello


greeting_jp = "こんにちは! 世界"  // Mixture of single byte string and multi byte strings.
str_len(greeting_jp)     // 8
str_subset(1,5)    // こんにちは
str_subset(8,9)    // 世界

```

### Regular expression support

Currently, the regular expression objects are not generated during run time but only in parsing phase. In other words, regular expressions can be generated only in regular expression literal. Regular expression syntax is re/pattern/ , which start with re/ and end with /. 

The regular expression engine used internally is Onigmo ( https://github.com/k-takata/Onigmo ). The source code must be written in UTF8, because regular expression object is programmed to assume UTF8 for its input (in libsailr).

As the following code shows, backreference (outside of the pattern) can be used via regexp_matched() function.

```
germany = re/(^Merc|^Porsche|^Volvo)/

if ( _rowname_ =~ germany ) {
    country = "Germany" ;
    type = rexp_matched(1);
}
```


### Date support

In Sailr(libsair), there are no specific types for dates, and dates are stored as just integers in Sailr. The number represents the day from Jan 1st 1970, which is the same as UNIX date. (Be careful that SAS software also store date information as numbers, but their starting date is 1/1/1960. https://support.sas.com/resources/papers/proceedings12/244-2012.pdf )

Internally for date calculation, date library is used ( "date.h" is a header-only library which builds upon <chrono>. https://github.com/HowardHinnant/date ).

```
date_ymd( year, month, day )
date_ym_weekday_nth( year, month, weekday, nth )  // weekday should be "Sun", "Mon" ... 
date_add_n_years( unix_date, years )
date_add_n_months( unix_date, months )
date_add_n_days( unix_date, days )
date_format( unix_date, str_format )  // This format should follow C++'s std::chrono::format. "%m/%d/%y" or "%Y-%m-%d" is popular.
```

### DataSailr specific functions

As already mentioned, DataSailr script is an extended version of Sailr script, and DataSailr has additional built-in functions to Sailr. Functions of Sailr script only conduct calculation within each row because that is what libsailr does. However, there are cases where users need to discard rows or output multiple rows from one row. DataSailr implements additional these functions using libsailr's external function mechanism. Currently, push!() and discard!() functions are provided using this mechanism. These functions can change the number of rows of input dataset, unlike libsailr's original built-in functions.


* push!() function

push!() function changes dataframe structure, for which reason it is not (cannot be) implemented as a libsailr internal function. What push!() function does is to push the current variable values onto the result dataframe. Even after push!() function is executed, the script execution continues on the same input row. This results in creating multiple rows from one input row.

As the following example shows, it is useful to convert wide format dataframe into long format.

```
wide_df = data.frame(
  subj = c("Tom", "Mary", "Jack"),
  t0 = c( 50, 42, 80),
  t1 = c( 48, 42, 75),
  t2 = c( 46, 44, 72),
  t3 = c( 42, 42, 73)
)

code = "
  OriRowNum = _n_
  subj = subj

  time = 0
  bw = t0
  push!()

  time = 1
  bw = t1 
  push!()

  time = 2
  bw = t2
  push!()

  time = 3
  bw = t3
"

wide_df

datasailr::sail(wide_df , code = code)
```


```
> wide_df

  subj t0 t1 t2 t3
1  Tom 50 48 46 42
2 Mary 42 42 44 42
3 Jack 80 75 72 73

> datasailr::sail(wide_df , code = code)

   subj OriRowNum time bw
1  Tom          1    0 50
2  Tom          1    1 48
3  Tom          1    2 46
4  Tom          1    3 42
5  Mary         2    0 42
6  Mary         2    1 42
7  Mary         2    2 44
8  Mary         2    3 42
9  Jack         3    0 80
10 Jack         3    1 75
11 Jack         3    2 72
12 Jack         3    3 73
```


* discard!() function

discard!() function drops the current row from the result dataframe. A common usage of this function is to filter out specific rows using if condition.

```
data(mtcars)

datasailr::sail(mtcars, code='
n = _n_
cyl = cyl
if(cyl != 4){
  discard!(); // Filter out rows that do not have 4 for cyl. i.e. Keep rows with 4 for cyl.
}
')
```


## For SAS software users

List of major differences between DataSailr and DATA Step in SAS software.

1. DataSailr if-else statements use curly brackets { }. SAS software uses if-then-else style.
2. DataSailr if-else curly brackets can hold multiple statements. SAS software needs to use do statement for multiple statements. i.e. if-then-do-end style. 
3. OUTPUT statement is not supported in DataSailr yet.
4. RETAIN statement is not supported and there is no plan to support it in the future.
5. The default character encoding of DataSailr is UTF8. 
6. Not only semicolon but also new line means the end of each statement in DataSailr. In SAS software, semicolons are always required.
7. (CAUTION) Date starts from Jan. 1st 1970, which is same as UNIX date. SAS software also stores date as number that is a day from 1/1/1960. Meaning NEVER DIRECTLY USE SAS NUMBER FOR DATES in DataSailr(libsailr).
8. Column name does not have label (In SAS software, dataset column can have a description called label.)



## For dplyr users

DataSailr does not compete with dplyr. DataSailr focuses on row by row manipulation, and dplyr has more functionalities about the whole dataframe manipulation (e.g. arrange(), select(), group_by() & summarize()).

DataSailr's special feature is that you can write multiple row level manipulations in one code, though dplyr enables it with a chain of functions.


|                        |     DataSailr                  |     dplyr            |
-------------------------|--------------------------------|------------------------ 
| How to manipulate data | Apply a single DataSailr code  | Apply multiple functions using (%>%) |
| Create new column      | assign value to new variable   | mutate()             |
| Keep some columns      | (Not implemented yet)          | select()             |
| Keep some rows         | discard!() drops rows          | filter()             |
| Summarize columns      | No use for this purpose        | summarize()          |
| Sort rows              | No use for this purpose        | arrange()            |
| Regular expression     | Built-in                       | Partially available with another R package |
| Available functions    | Can call only DataSailr functions  | Can call R functions |
| Convert wide to long format | push!() function          | (not dplyr but reshape2 package) |
| Convert long to wide format | (Not implemented yet)     | (not dplyr but reshape2 package) |



## For developers

### Architecture summary

The core engine of this package is libsailr. Libsailr takes Sailr script and does arithmetic calculations and string/character manipulations. This R package wraps libsair, and passes values of each row to the libsailr library. The results of each row will be returned as each row of new dataframe.

DataSailr internals are described in the official website.


## Community Guidelines

### When you wish to contribute to the software

This project's repository is on [Github](https://github.com/niceume/datasailr). Currently, for both development and releases, only master branch is used. 

Please feel free to fork the branch, and please send pull requests to me.


### When you report issues or problems with the software 

Please report issues or problems on Github.


### When you seek support

For support, there are online resources, a duscussion board and an email.

* You can read documents on [Official website](https://datasailr.io).
* You can ask questions on [Github discussion board](https://github.com/niceume/datasailr/discussions)
* You can send an email to the maintainer. (email address is described below)



## Contact

Your feedback is welcome. 

Maintainer: Toshi Umehara <toshi@niceume.com>


# HISTORY

## Version 0.8.9

* Bug fix in vec_list_push_cloned_row() function.

## Version 0.8.8

* This package now can be successfully built on FreeBSD.
* README is updated.
* Only sail() function is now exported. test_sail() function is now changed to an internal function.

## Version 0.8.7

* Thanksgiving Day Release (11/26/2020)
* In configure.ac, Onigmo's configure script is changed to be called with --disable-maintainer-mode, by which users do not re-run autoconf. (In some situations where timestamp information is not correctly preserved, autoconf tries to run again on user machines, which may cause aclocal missing error.)
    + To enable/disable maintainer mode, Onigmo's configure.ac needs to have 'AM_MAINTAINER_MODE([enable])'. DataSailr packaging toolchain modifies Onigmo's configure.ac to have this line.
* README and documents now use the term 'DataSailr script'.
* Options for configure script, --enable-libsailr-debug and --enable-datasailr-debug, are now available.
    + R CMD INSTALL can pass these options to configure with '--configure-args=' option.
    + By setting these options, -DDEBUG options are passed to compilers and DEBUG macro variable is defined.
* Fixes
    + When assignment does not happen in script, an error occured. However, there are cases where only discard!() and push!() are executed, and in such cases processed dataset should be returned if fullData = T. This is fixed.


## Version 0.8.6

* Update to follow libsailr's external function availability.
* push!() function is implemented as a DataSailr extended function.
    + It is implemented as an external function of libsailr, because it does more than just row level calculation. Libsailr cannot handle this. At DataSailr level, dataframe structure can be modified, and new rows can be added to result dataframe.
    + The function name, push!, comes from two reasons.
        + It is important to clarify that this function is different from Sailr built-in functions. This function works for more than each row level. It changes dataframe structures, and it has ! at the end.
        + This function does not stop but suspends libsailr vm execution, and adds or pushes rows. Therefore it is named as push.
    + The common usage of this function is converting data of wide format into long format.
* discard!() function is implemented as a DataSailr extended function.
    + discard!() function drops the current row from the result dataframe.
    + A common usage of this function is to filter out specific rows using if condition.


## Version 0.8.5

* Minor update
    + Warning with gcc link time optimization (on Fedora) is resolved.

## Version 0.8.4

* Cross-platform compilation is improved.
* The core processing library, libsailr has become able to report runtime errors with correspoding script location, and DataSailr supports it.
* Fixes
    + Variables that appear only on RHS of script crashed the program, which is fixed and the program stops safely.
    + String objects and regular expression objects were not freed in some cases. Memory leaks are fixed. 
    + NA in StringVector is treated as NULL, which results in passing an empty string to libsailr.


## Version 0.8.3

* Memory leaks fixed at various points that are pointed out after CRAN submission. (Apr.5-13 2020)
* Check whether parsing succeded or not. Function stops if the Sailr script has a wrong syntax. (Apr.14 2020)


## Version 0.8.2

* CRAN release candidate (Mar. 12 2020)
    + Biarch: true in DESCRIPTION.
        + In windows environment, if you have configure.win file, cross-compilation is not executed. This Biarch option solves it.
        + (ref) https://community.rstudio.com/t/configure-win-and-cran-submission/24684/4?fbclid=IwAR3RIMNzABwyBicz6hSPG6m1hmurtEzLGPjgMTFsjXJVkHEqBdW6wYFPd-Q
    + Package C++ main file (data_sailr_cpp_main.cpp) is updated for minor warnings.
        + Rcpp vectors' size() function returns singed int. Chage my code to be compatible with it.
    + Pass R defined makefile variables to Onigmo and libsailr makefiles
        + Makevars.win is updated.

* Minor updates (Mar. 15 2020)
    + Typos are corrected in README.md
    + Also, the following files are updated to resolve CRAN submission problems.
    + library() functions are removed.
        + These are used especially in test functions.
    + Title and description in DESCRIPTION are updated.

* Minor updates (Mar. 18 2020)
    + DESCRIPTION file is updated
    + printf() functions are enclosed within IF_DEBUG() macro in src/data_sailr_cpp_main.cpp
        + Including stdout/stderr related functions is not allowed in CRAN package.
    + bin directory is renamed to exec

* Minor updates (Mar. 22 2020)
    + DESCRIPTION file is updated. (Encoding and SystemRequirements)
    + Example script is added to package document.

* Factor class vector in dataframe is supported. They are now dealt as CharactorVector. (Apr. 1 2020)
    + In this implementation, factor levels is supposed to be CharacterVector.


## Version 0.8.1

* Beta release

## Version 0.8

* First commit
* Thanksgiving Day Release (11/28/2019)
* Christmas Update (12/25/2019) to release to the public
* Updating to submit to CRAN (01/11/2020)
    + The variables updated by this package have character types, not factors.
    + New argument, rowname, is added to sail() function.
    + libsailr is updated.
* Change following the libsailr API update. (01/19/2020)
* cleanup script is updated (2/4/2020)
* Resolve warnings for submitting CRAN (2/5/2020)
* Resolve windows compilation failures. (avoid autotools dependency during installation) (2/8/2020)
* configure.ac updates to cope with autotools nonexistent system. (2/9/2020)
* Makevars.win is updated. To build onigmo, win32/Makefile.mingw is now used. 
    + This leads to successfully build on rhub windows environment. (3/8/2020)


## Version 0.0 (Birth)

* Package skeleton was created. (11/15/2018)
   + The original name was RCppCalc, which was intended only for arithmetic calculation.

\name{adult_total}
\docType{data}
\alias{adult_total}
\title{Census income dataset}
\description{
"Census Income" dataset from UCI machine learning repository. 
}
\usage{data(adult_total)}
\format{
>50K, <=50K.
age: continuous.
workclass: Private, Self-emp-not-inc, Self-emp-inc, Federal-gov, Local-gov, State-gov, Without-pay, Never-worked.
fnlwgt: continuous.
education: Bachelors, Some-college, 11th, HS-grad, Prof-school, Assoc-acdm, Assoc-voc, 9th, 7th-8th, 12th, Masters, 1st-4th, 10th, Doctorate, 5th-6th, Preschool.
education-num: continuous.
marital-status: Married-civ-spouse, Divorced, Never-married, Separated, Widowed, Married-spouse-absent, Married-AF-spouse.
occupation: Tech-support, Craft-repair, Other-service, Sales, Exec-managerial, Prof-specialty, Handlers-cleaners, Machine-op-inspct, Adm-clerical, Farming-fishing, Transport-moving, Priv-house-serv, Protective-serv, Armed-Forces.
relationship: Wife, Own-child, Husband, Not-in-family, Other-relative, Unmarried.
race: White, Asian-Pac-Islander, Amer-Indian-Eskimo, Other, Black.
sex: Female, Male.
capital-gain: continuous.
capital-loss: continuous.
hours-per-week: continuous.
native-country: United-States, Cambodia, England, Puerto-Rico, Canada, Germany, Outlying-US(Guam-USVI-etc), India, Japan, Greece, South, China, Cuba, Iran, Honduras, Philippines, Italy, Poland, Jamaica, Vietnam, Mexico, Portugal, Ireland, France, Dominican-Republic, Laos, Ecuador, Taiwan, Haiti, Columbia, Hungary, Guatemala, Nicaragua, Scotland, Thailand, Yugoslavia, El-Salvador, Trinadad&Tobago, Peru, Hong, Holand-Netherlands.
}
\source{
  https://archive.ics.uci.edu/ml/datasets/Adult
}
\references{
  https://archive.ics.uci.edu/ml/datasets/Adult
}
\keyword{datasets}
%File sail.Rd

\name{sail}
\alias{sail}
\title{Run 'DataSailr' script and manipulate dataframe }
\description{
  datasailr::sail() is the main function of this 'DataSailr' package. This function takes a data.frame as the 1st argument and a 'DataSailr' script as the 2nd argument. 
}
\usage{
  sail( df, code, fullData = TRUE, rowname = "_rowname_", stringsAsFactors = FALSE)
}
\arguments{
  \item{df}{ data.frame to be processed.}
  \item{code}{ 'DataSailr' script that instructs how to manipulate values and strings.}
  \item{fullData}{ When fullData is set FALSE, only the columns corresponding to left-hand-side(LHS) variables in 'DataSailr' scripts are returned. (e.g. In 'DataSailr' script, bmi = weight / height / height , the variable bmi is the LHS.) If set TRUE, all the columns are returned. }
  \item{rowname}{ Dataframe's row names are accessed via the variable name specified by this option. The default variable name is _rowname_, meaning that _rowname_ represents dataframe's row names in 'DataSailr' script, though this does not allow assignment. Even if you assign some value to this variable, they are not reflected in the row names of the resultant dataframe. }
  \item{stringsAsFactors}{ Columns of character vectors (= string vectors) that appear as left-hand-side(LHS) variables in 'DataSailr' scripts can be returned as either character vectors or factor vectors. If stringsAsFactors is set to FALSE(default), character columns are returned as character vectors. If stringsAsFactors is set to TRUE, character columns are returned as factor vectors. }
}
\value{
  If fullData argument is set FALSE, return value is a data.frame. If fullData is TRUE, the result is created by cbind(), which concatenates the original input and newly created dataframe. For example, if the original input is 'data.table', the return value is 'data.table'.
}
\examples{
library(datasailr)
data(iris)

iris_code = '
iris_type = .
if ( Species == "setosa" ) { 
    type = 1
    type_char = "A"
}else if( Species == "versicolor" ) { 
    type = 2
    type_char = "B"
}else if(Species == "virginica" ) { 
    type = 3
    type_char = "C"
}
'
sail(iris, iris_code)
}

%File datasailr-package.Rd

\name{datasailr-package}
\alias{datasailr-package}
\docType{package}
\title{
  Row by Row Data Processing Tool, Using 'DataSailr' Script
}
\description{
  A row by row data processing tool. You can write data manipulation code in 'DataSailr' script which is specially intended for the data manipulation. The package uses 'libsailr' (C/C++ library) for its 'DataSailr' code parsing and its execution.
}
\details{
  The main function, datasailr::sail(), takes dataframe as its 1st argument and 'DataSailr' script as its 2nd argument, and returns the processing results as dataframe. This package works with 'libsailr' library, which conducts arithmetic calculations and string/character manipulations following 'Sailr' script (which is a base version of 'DataSailr' script). 'DataSailr' package passes both 'DataSailr' script and values of each row of dataframe to the 'libsailr', and 'libsailr' processes those values of each row. The processing results are stored and are returned as dataframe.
}
\author{
Toshi Umehara, toshi@niceume.com.
Maintainer: Toshi Umehara <toshi@niceume.com>
}
\keyword{ package }
\seealso{
  \code{\link{sail}}.
}
\examples{
library(datasailr)
data(mtcars)

code_example = 'powerful = .
  if( 145 >= hp && hp > 0 ){
    powerful = 0
  }else if( hp > 145 ){
    powerful = 1
  }else{
    powerful = 2
  }

  germany = re/(^Merc|^Porsche|^Volvo)/
  japan = re/(^Mazda|^Honda|^Toyota)/

  if ( _rowname_ =~  germany ) { 
      country = "Germany" 
  }else if( _rowname_ =~ japan ){
      country = "Japan"
  }else {
      country = "Other"
  }
  company = rexp_matched(1)
'

sail(mtcars, code_example)
}

