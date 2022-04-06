# OmicsView

Demo System: http://omicsview.org/

OmicsView Supplementary: https://github.com/interactivereport/OmicsView/blob/main/OmicsView_Supp_v5.0.pdf

Installation Guide: http://omicsview.org/docs/

OmicsView: omics data analysis through interactive visual analytics

Fergal Casey<sup>1</sup>, Soumya Negi<sup>1</sup>, Joost Groot<sup>1</sup>, Jing Zhu<sup>1</sup>, Maria Zavodszky<sup>1</sup>, Derrick Cheng<sup>2</sup>, Dongdong Lin<sup>1</sup>, Sally John<sup>1</sup>, Michelle A. Penny<sup>1</sup>, David Sexton<sup>1</sup>, Baohong Zhang<sup>1</sup>

<sup>1</sup>Translational Biology, Research Development, Biogen, Inc., Cambridge, MA, 02142, USA  
<sup>2</sup>BioInfoRx, Inc., 510 Charmany Dr, Suite 275A, Madison, WI 53719, USA



Please contact Baohong Zhang (baohongz@gmail.com) if you have any question about this project.
# jstree

[jsTree](http://www.jstree.com/) is jquery plugin, that provides interactive trees. It is absolutely free, [open source](https://github.com/vakata/jstree) and distributed under the MIT license.

jsTree is easily extendable, themable and configurable, it supports HTML & JSON data sources, AJAX & async callback loading.

jsTree functions properly in either box-model (content-box or border-box), can be loaded as an AMD module, and has a built in mobile theme for responsive design, that can easily be customized. It uses jQuery's event system, so binding callbacks on various events in the tree is familiar and easy.

You also get:
 * drag & drop support
 * keyboard navigation
 * inline edit, create and delete
 * tri-state checkboxes
 * fuzzy searching
 * customizable node types

_Aside from this readme you can find a lot more info on [jstree.com](http://www.jstree.com) & [the discussion group](https://groups.google.com/forum/#!forum/jstree)_.

---

<!-- MarkdownTOC depth=0 autolink=true bracket=round -->

- [Getting Started](#getting-started)
  - [Include all neccessary files](#include-all-neccessary-files)
  - [Populating a tree using HTML](#populating-a-tree-using-html)
  - [Populating a tree using an array \(or JSON\)](#populating-a-tree-using-an-array-or-json)
    - [The required JSON format](#the-required-json-format)
  - [Populating the tree using AJAX](#populating-the-tree-using-ajax)
  - [Populating the tree using AJAX and lazy loading nodes](#populating-the-tree-using-ajax-and-lazy-loading-nodes)
  - [Populating the tree using a callback function](#populating-the-tree-using-a-callback-function)
- [Working with events](#working-with-events)
- [Interacting with the tree using the API](#interacting-with-the-tree-using-the-api)
- [More on configuration](#more-on-configuration)
- [Plugins](#plugins)
  - [checkbox](#checkbox)
  - [contextmenu](#contextmenu)
  - [dnd](#dnd)
  - [massload](#massload)
  - [search](#search)
  - [sort](#sort)
  - [state](#state)
  - [types](#types)
  - [unique](#unique)
  - [wholerow](#wholerow)
  - [More plugins](#more-plugins)
- [License & Contributing](#license--contributing)

<!-- /MarkdownTOC -->


---

## Getting Started

### Include all neccessary files
To get started you need 3 things in your page:
 1. jQuery (anything above 1.9.1 will work)
 2. A jstree theme (there is only one theme supplied by default)
 3. The jstree source file

```html
<script src="//cdnjs.cloudflare.com/ajax/libs/jquery/3.1.1/jquery.min.js"></script>

<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/themes/default/style.min.css" />
<script src="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/jstree.min.js"></script>
```

_If you decide to host jstree yourself - the files are located in the `dist` folder. You can safely ignore the `dist/libs` folder._

---

### Populating a tree using HTML

Now we are all set to create a tree, inline HTML is the easiest option (suitable for menus). All you need to do is select a node (using a jQuery selector) and invoke the `.jstree()` function to let jstree know you want to render a tree inside the selected node. `$.jstree.create(element)` can be used too.

```html
<div id="container">
  <ul>
    <li>Root node
      <ul>
        <li>Child node 1</li>
        <li>Child node 2</li>
      </ul>
    </li>
  </ul>
</div>
<script>
$(function() {
  $('#container').jstree();
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/)

_You can add a few options when rendering a node using a data-attribute (note the quotes):_
```html
<li data-jstree='{ "selected" : true, "opened" : true }'>Root node ...
```

---

### Populating a tree using an array (or JSON)

Building trees from HTML is easy, but it is not very flexible, inline JS data is a better option:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        { "text" : "Root node", "children" : [
            { "text" : "Child node 1" },
            { "text" : "Child node 2" }
          ]
        }
      ]
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4478/)

Unlike the previous simple HTML example, this time the `.jstree()` function accepts a config object.

For now it is important to note that jstree will try to parse any data you specify in the  `core.data` key and use it to create a tree. As seen in the previous example, if this key is missing jstree will try to parse the inline HTML of the container.

#### The required JSON format

The data you use must be in a specific format, each branch of the tree is represented by an object, which must at least have a `text` key. The `children` key can be used to add children to the branch, it should be an array of objects.

_Keep in mind, you can use a simple string instead of an object if all you need is node with the given text, the above data can be written as:_

```js
[ { "text" : "Root node", "children" : [ "Child node 1", "Child node 2" ] } ]
```

There are other available options for each node, only set them if you need them like:

 * `id` - makes if possible to identify a node later (will also be used as a DOM ID of the `LI` node). _Make sure you do not repeat the same ID in a tree instance (that would defeat its purpose of being a unique identifier and may cause problems for jstree)_.
 * `icon` - a string which will be used for the node's icon - this can either be a path to a file, or a className (or list of classNames), which you can style in your CSS (font icons also work).
 * `data` - this can be anything you want - it is metadata you want attached to the node - you will be able to access and modify it any time later - it has no effect on the visuals of the node.
 * `state` - an object specifyng a few options about the node:
   - `selected` - if the node should be initially selected
   - `opened` - if the node should be initially opened
   - `disabled` - if the node should be disabled
   - `checked` - __checkbox plugin specific__ - if the node should be checked (only used when `tie_selection` is `false`, which you should only do if you really know what you are doing)
   - `undetermined` - __checkbox plugin specific__ - if the node should be rendered in undetermined state (only used with lazy loading and when the node is not yet loaded, otherwise this state is automatically calculated).
 * `type` - __types plugin specific__ - the type of the nodes (should be defined in the types config), if not set `"default"` is assumed.
 * `li_attr` - object of values which will be used to add HTML attributes on the resulting `LI` DOM node.
 * `a_attr` - object of values which will be used to add HTML attributes on the resulting `A` node.

Here is a new demo with some of those properties set:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
          {
              "text" : "Root node",
              "state" : {"opened" : true },
              "children" : [
                  {
                    "text" : "Child node 1",
                    "state" : { "selected" : true },
                    "icon" : "glyphicon glyphicon-flash"
                  },
                  { "text" : "Child node 2", "state" : { "disabled" : true } }
              ]
        }
      ]
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4479/)

---

### Populating the tree using AJAX

Building off of the previous example, let's see how to have jstree make AJAX requests for you.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : {
        "url" : "//www.jstree.com/fiddle/",
        "dataType" : "json" // needed only if you do not supply JSON headers
      }
    }
  });
});
</script>
```

The server response is:
```json
[{
  "id":1,"text":"Root node","children":[
    {"id":2,"text":"Child node 1"},
    {"id":3,"text":"Child node 2"}
  ]
}]
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4480/)

Instead of a JS array, you can set `core.data` to a [jQuery AJAX config](http://api.jquery.com/jQuery.ajax/). 
jsTree will hit that URL, and provided you return properly formatted JSON it will be displayed.

_If you cannot provide proper JSON headers, set `core.data.dataType` to `"json"`._

The ids in the server response make it possible to identify nodes later (which we will see in the next few demos), but they are not required.

__WHEN USING IDS MAKE SURE THEY ARE UNIQUE INSIDE A PARTICULAR TREE__

---

### Populating the tree using AJAX and lazy loading nodes

Lazy loading means nodes will be loaded when they are needed. Imagine you have a huge amount of nodes you want to show, but loading them with a single request is way too much traffic. Lazy loading makes it possible to load nodes on the fly - jstree will perform AJAX requests as the user browses the tree.

Here we take our previous example, and lazy load the "Child node 1" node.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : {
        "url" : "//www.jstree.com/fiddle/?lazy",
        "data" : function (node) {
          return { "id" : node.id };
        }
      }
    }
  });
});
</script>
```

The initial server response is:
```json
[{
  "id":1,"text":"Root node","children":[
    {"id":2,"text":"Child node 1","children":true},
    {"id":3,"text":"Child node 2"}
  ]
}]
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4481/)

Now to focus on what is different. First off the `"data"` config option of the data object. If you check with jQuery, it is supposed to be a string or an object. But jstree makes it possible to set a function.

Each time jstree needs to make an AJAX call this function will be called and will receive a single parameter - the node that is being loaded. The return value of this function will be used as the actual `"data"` of the AJAX call. To understand better open up the demo and see the requests go off in the console.

You will notice that the first request goes off to:
`http://www.jstree.com/fiddle?lazy&id=#`
`#` is the special ID that the function receives when jstree needs to load the root nodes.

Now go ahead and open the root node - two children will be shown, but no request will be made - that is because we loaded those children along with the first request.

Onto the next difference - "Child node 1" appears closed - that is because in the data we supplied `true` as the `"children"` property of this node (you can see it in the server response). This special value indicated to jstree, that it has to lazy load the "Child node 1" node.

Proceed and open this node - you will see a next request fire off to:
`http://www.jstree.com/fiddle?lazy&id=2`
ID is set to `2` because the node being loaded has an ID of `2`, and we have configured jstree to send the node ID along with the AJAX request (the `data` function).

The server response is:
```json
["Child node 3","Child node 4"]
```

_You can also set `"url"` to a function and it works exactly as with `"data"` - each time a request has to be made, jstree will invoke your function and the request will go off to whatever you return in this function. This is useful when dealing with URLs like: `http://example.com/get_children/1`._

### Populating the tree using a callback function

Sometimes you may not want jsTree to make AJAX calls for you - you might want to make them yourself, or use some other method of puplating the tree. In that case you can use a callback function.

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : function (node, cb) {
        if(node.id === "#") {
          cb([{"text" : "Root", "id" : "1", "children" : true}]);
        }
        else {
          cb(["Child"]);
        }
      }
    }
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4482/)

As you can see your function will receive two arguments - the node whose children need to be loaded and a callback function to call with the data once you have it. The data follows the same familiar JSON format and lazy loading works just as with AJAX (as you can see in the above example).

---

## Working with events

jstree provides a lot of events to let you know something happened with the tree. The events are the same regardless of how you populate the tree.
Let's use the most basic event `changed` - it fires when selection on the tree changes:

```html
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        {"id" : 1, "text" : "Node 1"},
        {"id" : 2, "text" : "Node 2"},
      ]
    }
  });
  $('#container').on("changed.jstree", function (e, data) {
    console.log("The selected nodes are:");
    console.log(data.selected);
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4483/)

All jstree events fire in a special `".jstree"` namespace - this is why we listen for `"changed.jstree"`. The handler itself receives one additional parameter - it will be populated with all you need to know about the event that happened. In this case `data.selected` is an array of selected node IDs (please note, that if you have not specified IDs they will be autogenerated).

Let's extend this a bit and log out the text of the node instead of the ID.

```js
$('#container').on("changed.jstree", function (e, data) {
  console.log(data.instance.get_selected(true)[0].text);
  console.log(data.instance.get_node(data.selected[0]).text);
});
```

The two rows above achieve exactly the same thing - get the text of the first selected node.

In the `data` argument object you will always get an `instance` key - that is a reference to the tree instance, so that you can easily invoke methods.

__All available functions and events are documented in the API docs__

---

## Interacting with the tree using the API

We scratched the surface on interacting with the tree in the previous example. Let's move on to obtaining an instance and calling a method on this instance:

```html
<button>Select node 1</button>
<div id="container"></div>
<script>
$(function() {
  $('#container').jstree({
    'core' : {
      'data' : [
        {"id" : 1, "text" : "Node 1"},
        {"id" : 2, "text" : "Node 2"},
      ]
    }
  });
  $('button').on("click", function () {
    var instance = $('#container').jstree(true);
    instance.deselect_all();
    instance.select_node('1');
  });
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4484/)

The above example shows how to obtain a reference to a jstree instance (again with a selector, but this time instead of a config, we pass a boolean `true`), and call a couple of methods - the latter one is selecting a node by its ID.

Methods can also be invoked like this:

```js
$('#container').jstree("select_node", "1");
```

__All available functions and events are documented in the API docs__

## More on configuration

We already covered the config object in general (when we specified inline & AJAX data sources).

```js
$("#tree").jstree({ /* config object goes here */ });
```

Each key in the config object corresponds to a plugin, and the value of that key is the configuration for that plugin. There are also two special keys `"core"` and `"plugins"`:
 * `"core"` stores the core configuration options
 * `"plugins"` is an array of plugin names (strings) you want active on the instance

When configuring you only need to set values that you want to be different from the defaults.

__All config options and defaults are documented in the API docs__

```js
$("#tree").jstree({
  "core" : { // core options go here
    "multiple" : false, // no multiselection
    "themes" : {
      "dots" : false // no connecting dots between dots
    }
  },
  "plugins" : ["state"] // activate the state plugin on this instance
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4485/)

We will cover all plugins further down.

__Keep in mind by default all modifications to the structure are prevented - that means drag'n'drop, create, rename, delete will not work unless you enable them.__

```js
$("#tree").jstree({
  "core" : {
    "check_callback" : true, // enable all modifications
  },
  "plugins" : ["dnd","contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4486/)

`"core.check_callback"` can also be set to a function, that will be invoked every time a modification is about to happen (or when jstree needs to check if a modification is possible). If you return `true` the operation will be allowed, a value of `false` means it will not be allowed. The possible operation you can expect are `create_node`, `rename_node`, `delete_node`, `move_node` and `copy_node`. The `more` parameter will contain various information provided by the plugin that is invoking the check. For example the DND plugin will provide an object containing information about the move or copy operation that is being checked - is it a multi tree operation, which node is currently hovered, where the insert arrow is pointing - before, after or inside, etc.

```js
$("#tree").jstree({
  "core" : {
    "check_callback" : function (operation, node, parent, position, more) {
      if(operation === "copy_node" || operation === "move_node") {
        if(parent.id === "#") {
          return false; // prevent moving a child above or below the root
        }
      },
      return true; // allow everything else
    }
  },
  "plugins" : ["dnd","contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4487/)

The `more` parameter you receive contains other information related to the check being performed.

__For example__: `move_node` & `copy_node` checks will fire repeatedly while the user drags a node, if the check was triggered by the `dnd` plugin `more` will contain a `dnd` key, which will be set to `true`.
You can check for `more.dnd` and only perform a certain action if `dnd` triggered the check.
If you only want to perform an operation when a node is really about to be dropped check for `more.core`.

## Plugins

jsTree comes with a few plugin bundled, but they will only modify your tree if you activate them using the `"plugins"` config option. Here is a brief description of each plugin. You can read more on the available config options for each plugin in the API docs.

### checkbox
Renders a checkbox icon in front of each node, making multiselection easy. It also has a "tri-state" option, meaning a node with some of its children checked will get a "square" icon.

_Keep in mind that if any sort of cascade is enabled, disabled nodes may be checked too (not by themselves, but for example when a parent of a disabled node is checked and selection is configured to cascade down)._

```js
$("#tree").jstree({
  "plugins" : ["checkbox"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4488/)

### contextmenu
Makes it possible to right click nodes and shows a list of configurable actions in a menu.

```js
$("#tree").jstree({
  "core" : { "check_callback" : true }, // so that modifying operations work
  "plugins" : ["contextmenu"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4489/)

### dnd
Makes it possible to drag and drop tree nodes and rearrange the tree.

```js
$("#tree").jstree({
  "core" : { "check_callback" : true }, // so that operations work
  "plugins" : ["dnd"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4490/)

### massload
Makes it possible to load multiple nodes in a single go (for a lazy loaded tree).

```js
$("#tree").jstree({
  "core" : {
    "data" : { .. AJAX config .. }
  },
  "massload" : {
    "url" : "/some/path",
    "data" : function (nodes) {
      return { "ids" : nodes.join(",") };
    }
  },
  "plugins" : [ "massload", "state" ]
});
```

### search
Adds the possibility to search for items in the tree and show only matching nodes. It also has AJAX / callback hooks, so that search will work on lazy loaded trees too.

```html
<form id="s">
  <input type="search" id="q" />
  <button type="submit">Search</button>
</form>
<script>
$("#container").jstree({
  "plugins" : ["search"]
});
$("#s").submit(function(e) {
  e.preventDefault();
  $("#container").jstree(true).search($("#q").val());
});
</script>
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4491/)

### sort
Automatically arranges all sibling nodes according to a comparison config option function, which defaults to alphabetical order.

```js
$("#tree").jstree({
  "plugins" : ["sort"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4492/)

### state
Saves all opened and selected nodes in the user's browser, so when returning to the same tree the previous state will be restored.

```js
$("#tree").jstree({
  // the key is important if you have multiple trees in the same domain
  "state" : { "key" : "state_demo" },
  "plugins" : ["state"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4493/)

### types
Makes it possible to add a "type" for a node, which means to easily control nesting rules and icon for groups of nodes instead of individually. To set a node type add a type property to the node structure.

```js
$("#tree").jstree({
  "types" : {
    "default" : {
      "icon" : "glyphicon glyphicon-flash"
    },
    "demo" : {
      "icon" : "glyphicon glyphicon-ok"
    }
  },
  "plugins" : ["types"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4494/)

### unique
Enforces that no nodes with the same name can coexist as siblings - prevents renaming and moving nodes to a parent, which already contains a node with the same name.

```js
$("#tree").jstree({
  "plugins" : ["unique"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4495/)

### wholerow
Makes each node appear block level which makes selection easier. May cause slow down for large trees in old browsers.

```js
$("#tree").jstree({
  "plugins" : ["wholerow"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4496/)

### More plugins
If you create your own plugin (or download a 3rd party one) you must include its source on the page and list its name in the `"plugins"` config array.

```js
// conditional select
(function ($, undefined) {
  "use strict";
  $.jstree.defaults.conditionalselect = function () { return true; };
  $.jstree.plugins.conditionalselect = function (options, parent) {
    this.activate_node = function (obj, e) {
      if(this.settings.conditionalselect.call(this, this.get_node(obj))) {
        parent.activate_node.call(this, obj, e);
      }
    };
  };
})(jQuery);
$("#tree").jstree({
  "conditionalselect" : function (node) {
    return node.text === "Root node" ? false : true;
  },
  "plugins" : ["conditionalselect"]
});
```

[view result](http://jsfiddle.net/vakata/2kwkh2uL/4497/)

As seen here when creating a plugin you can define a default config, add your own functions to jstree, or override existing ones while maintaining the ability to call the overridden function.

## License & Contributing

_Please do NOT edit files in the "dist" subdirectory as they are generated via grunt. You'll find source code in the "src" subdirectory!_

If you want to you can always [donate a small amount][paypal] to help the development of jstree.
[paypal]: https://www.paypal.com/cgi-bin/webscr?cmd=_xclick&business=paypal@vakata.com&currency_code=USD&amount=&return=http://jstree.com/donation&item_name=Buy+me+a+coffee+for+jsTree

Copyright (c) 2014 Ivan Bozhanov (http://vakata.com)

Licensed under the [MIT license](http://www.opensource.org/licenses/mit-license.php).
This documentation is based on our [OAI specification](https://github.com/sendgrid/sendgrid-oai).

# INITIALIZATION

```php
// If you are using Composer
require 'vendor/autoload.php';


$apiKey = getenv('SENDGRID_API_KEY');
$sg = new \SendGrid($apiKey);
```

# Table of Contents

* [ACCESS SETTINGS](#access_settings)
* [ALERTS](#alerts)
* [API KEYS](#api_keys)
* [ASM](#asm)
* [BROWSERS](#browsers)
* [CAMPAIGNS](#campaigns)
* [CATEGORIES](#categories)
* [CLIENTS](#clients)
* [CONTACTDB](#contactdb)
* [DEVICES](#devices)
* [GEO](#geo)
* [IPS](#ips)
* [MAIL](#mail)
* [MAIL SETTINGS](#mail_settings)
* [MAILBOX PROVIDERS](#mailbox_providers)
* [PARTNER SETTINGS](#partner_settings)
* [SCOPES](#scopes)
* [SENDERS](#senders)
* [STATS](#stats)
* [SUBUSERS](#subusers)
* [SUPPRESSION](#suppression)
* [TEMPLATES](#templates)
* [TRACKING SETTINGS](#tracking_settings)
* [USER](#user)
* [WHITELABEL](#whitelabel)


<a name="access_settings"></a>
# ACCESS SETTINGS

## Retrieve all recent access attempts

**This endpoint allows you to retrieve a list of all of the IP addresses that recently attempted to access your account either through the User Interface or the API.**

IP Access Management allows you to control which IP addresses can be used to access your account, either through the User Interface or the API. There is no limit to the number of IP addresses that you can add to your whitelist. It is possible to remove your own IP address from the whitelist, thus preventing yourself from accessing your account.

For more information, please see our [User Guide](http://sendgrid.com/docs/User_Guide/Settings/ip_access_management.html).

### GET /access_settings/activity


```php
$query_params = json_decode('{"limit": 1}');
$response = $sg->client->access_settings()->activity()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add one or more IPs to the whitelist

**This endpoint allows you to add one or more IP addresses to your IP whitelist.**

When adding an IP to your whitelist, include the IP address in an array. You can whitelist one IP at a time, or you can whitelist multiple IPs at once.

IP Access Management allows you to control which IP addresses can be used to access your account, either through the User Interface or the API. There is no limit to the number of IP addresses that you can add to your whitelist. It is possible to remove your own IP address from the whitelist, thus preventing yourself from accessing your account.

For more information, please see our [User Guide](http://sendgrid.com/docs/User_Guide/Settings/ip_access_management.html).

### POST /access_settings/whitelist


```php
$request_body = json_decode('{
  "ips": [
    {
      "ip": "192.168.1.1"
    },
    {
      "ip": "192.*.*.*"
    },
    {
      "ip": "192.168.1.3/32"
    }
  ]
}');
$response = $sg->client->access_settings()->whitelist()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a list of currently whitelisted IPs

**This endpoint allows you to retrieve a list of IP addresses that are currently whitelisted.**

IP Access Management allows you to control which IP addresses can be used to access your account, either through the User Interface or the API. There is no limit to the number of IP addresses that you can add to your whitelist. It is possible to remove your own IP address from the whitelist, thus preventing yourself from accessing your account.

For more information, please see our [User Guide](http://sendgrid.com/docs/User_Guide/Settings/ip_access_management.html).

### GET /access_settings/whitelist


```php
$response = $sg->client->access_settings()->whitelist()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Remove one or more IPs from the whitelist

**This endpoint allows you to remove one or more IPs from your IP whitelist.**

You can remove one IP at a time, or you can remove multiple IP addresses.

IP Access Management allows you to control which IP addresses can be used to access your account, either through the User Interface or the API. There is no limit to the number of IP addresses that you can add to your whitelist. It is possible to remove your own IP address from the whitelist, thus preventing yourself from accessing your account.

For more information, please see our [User Guide](http://sendgrid.com/docs/User_Guide/Settings/ip_access_management.html).

### DELETE /access_settings/whitelist


```php
$request_body = json_decode('{
  "ids": [
    1,
    2,
    3
  ]
}');
$response = $sg->client->access_settings()->whitelist()->delete($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a specific whitelisted IP

**This endpoint allows you to retreive a specific IP address that has been whitelisted.**

You must include the ID for the specific IP address you want to retrieve in your call.

IP Access Management allows you to control which IP addresses can be used to access your account, either through the User Interface or the API. There is no limit to the number of IP addresses that you can add to your whitelist. It is possible to remove your own IP address from the whitelist, thus preventing yourself from accessing your account.

For more information, please see our [User Guide](http://sendgrid.com/docs/User_Guide/Settings/ip_access_management.html).

### GET /access_settings/whitelist/{rule_id}


```php
$rule_id = "test_url_param";
$response = $sg->client->access_settings()->whitelist()->_($rule_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Remove a specific IP from the whitelist

**This endpoint allows you to remove a specific IP address from your IP whitelist.**

When removing a specific IP address from your whitelist, you must include the ID in your call.

IP Access Management allows you to control which IP addresses can be used to access your account, either through the User Interface or the API. There is no limit to the number of IP addresses that you can add to your whitelist. It is possible to remove your own IP address from the whitelist, thus preventing yourself from accessing your account.

For more information, please see our [User Guide](http://sendgrid.com/docs/User_Guide/Settings/ip_access_management.html).

### DELETE /access_settings/whitelist/{rule_id}


```php
$rule_id = "test_url_param";
$response = $sg->client->access_settings()->whitelist()->_($rule_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="alerts"></a>
# ALERTS

## Create a new Alert

**This endpoint allows you to create a new alert.**

Alerts allow you to specify an email address to receive notifications regarding your email usage or statistics.
* Usage alerts allow you to set the threshold at which an alert will be sent.
* Stats notifications allow you to set how frequently you would like to receive email statistics reports. For example, "daily", "weekly", or "monthly".

For more information about alerts, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/alerts.html).

### POST /alerts


```php
$request_body = json_decode('{
  "email_to": "example@example.com",
  "frequency": "daily",
  "type": "stats_notification"
}');
$response = $sg->client->alerts()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all alerts

**This endpoint allows you to retieve all of your alerts.**

Alerts allow you to specify an email address to receive notifications regarding your email usage or statistics.
* Usage alerts allow you to set the threshold at which an alert will be sent.
* Stats notifications allow you to set how frequently you would like to receive email statistics reports. For example, "daily", "weekly", or "monthly".

For more information about alerts, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/alerts.html).

### GET /alerts


```php
$response = $sg->client->alerts()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update an alert

**This endpoint allows you to update an alert.**

Alerts allow you to specify an email address to receive notifications regarding your email usage or statistics.
* Usage alerts allow you to set the threshold at which an alert will be sent.
* Stats notifications allow you to set how frequently you would like to receive email statistics reports. For example, "daily", "weekly", or "monthly".

For more information about alerts, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/alerts.html).

### PATCH /alerts/{alert_id}


```php
$request_body = json_decode('{
  "email_to": "example@example.com"
}');
$alert_id = "test_url_param";
$response = $sg->client->alerts()->_($alert_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a specific alert

**This endpoint allows you to retrieve a specific alert.**

Alerts allow you to specify an email address to receive notifications regarding your email usage or statistics.
* Usage alerts allow you to set the threshold at which an alert will be sent.
* Stats notifications allow you to set how frequently you would like to receive email statistics reports. For example, "daily", "weekly", or "monthly".

For more information about alerts, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/alerts.html).

### GET /alerts/{alert_id}


```php
$alert_id = "test_url_param";
$response = $sg->client->alerts()->_($alert_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete an alert

**This endpoint allows you to delete an alert.**

Alerts allow you to specify an email address to receive notifications regarding your email usage or statistics.
* Usage alerts allow you to set the threshold at which an alert will be sent.
* Stats notifications allow you to set how frequently you would like to receive email statistics reports. For example, "daily", "weekly", or "monthly".

For more information about alerts, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/alerts.html).

### DELETE /alerts/{alert_id}


```php
$alert_id = "test_url_param";
$response = $sg->client->alerts()->_($alert_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="api_keys"></a>
# API KEYS

## Create API keys

**This enpoint allows you to create a new random API Key for the user.**

A JSON request body containing a "name" property is required. If number of maximum keys is reached, HTTP 403 will be returned.

There is a limit of 100 API Keys on your account.

The API Keys feature allows customers to be able to generate an API Key credential which can be used for authentication with the SendGrid v3 Web API or the [Mail API Endpoint](https://sendgrid.com/docs/API_Reference/Web_API/mail.html).

See the [API Key Permissions List](https://sendgrid.com/docs/API_Reference/Web_API_v3/API_Keys/api_key_permissions_list.html) for a list of all available scopes.

### POST /api_keys


```php
$request_body = json_decode('{
  "name": "My API Key",
  "sample": "data",
  "scopes": [
    "mail.send",
    "alerts.create",
    "alerts.read"
  ]
}');
$response = $sg->client->api_keys()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all API Keys belonging to the authenticated user

**This endpoint allows you to retrieve all API Keys that belong to the authenticated user.**

The API Keys feature allows customers to be able to generate an API Key credential which can be used for authentication with the SendGrid v3 Web API or the [Mail API Endpoint](https://sendgrid.com/docs/API_Reference/Web_API/mail.html).

### GET /api_keys


```php
$query_params = json_decode('{"limit": 1}');
$response = $sg->client->api_keys()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update the name & scopes of an API Key

**This endpoint allows you to update the name and scopes of a given API key.**

A JSON request body with a "name" property is required.
Most provide the list of all the scopes an api key should have.

The API Keys feature allows customers to be able to generate an API Key credential which can be used for authentication with the SendGrid v3 Web API or the [Mail API Endpoint](https://sendgrid.com/docs/API_Reference/Web_API/mail.html).


### PUT /api_keys/{api_key_id}


```php
$request_body = json_decode('{
  "name": "A New Hope",
  "scopes": [
    "user.profile.read",
    "user.profile.update"
  ]
}');
$api_key_id = "test_url_param";
$response = $sg->client->api_keys()->_($api_key_id)->put($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update API keys

**This endpoint allows you to update the name of an existing API Key.**

A JSON request body with a "name" property is required.

The API Keys feature allows customers to be able to generate an API Key credential which can be used for authentication with the SendGrid v3 Web API or the [Mail API Endpoint](https://sendgrid.com/docs/API_Reference/Web_API/mail.html).

## URI Parameters

| URI Parameter   | Type  | Required?  | Description  |
|---|---|---|---|
|api_key_id |string | required | The ID of the API Key you are updating.|

### PATCH /api_keys/{api_key_id}


```php
$request_body = json_decode('{
  "name": "A New Hope"
}');
$api_key_id = "test_url_param";
$response = $sg->client->api_keys()->_($api_key_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve an existing API Key

**This endpoint allows you to retrieve a single api key.**

If the API Key ID does not exist an HTTP 404 will be returned.

### GET /api_keys/{api_key_id}


```php
$api_key_id = "test_url_param";
$response = $sg->client->api_keys()->_($api_key_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete API keys

**This endpoint allows you to revoke an existing API Key**

Authentications using this API Key will fail after this request is made, with some small propogation delay.If the API Key ID does not exist an HTTP 404 will be returned.

The API Keys feature allows customers to be able to generate an API Key credential which can be used for authentication with the SendGrid v3 Web API or the [Mail API Endpoint](https://sendgrid.com/docs/API_Reference/Web_API/mail.html).

## URI Parameters

| URI Parameter   | Type  | Required?  | Description  |
|---|---|---|---|
|api_key_id |string | required | The ID of the API Key you are deleting.|

### DELETE /api_keys/{api_key_id}


```php
$api_key_id = "test_url_param";
$response = $sg->client->api_keys()->_($api_key_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="asm"></a>
# ASM

## Create a new suppression group

**This endpoint allows you to create a new suppression group.**

Suppression groups, or unsubscribe groups, are specific types or categories of email that you would like your recipients to be able to unsubscribe from. For example: Daily Newsletters, Invoices, System Alerts.

The **name** and **description** of the unsubscribe group will be visible by recipients when they are managing their subscriptions.

Each user can create up to 25 different suppression groups.

### POST /asm/groups


```php
$request_body = json_decode('{
  "description": "Suggestions for products our users might like.",
  "is_default": true,
  "name": "Product Suggestions"
}');
$response = $sg->client->asm()->groups()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve information about multiple suppression groups

**This endpoint allows you to retrieve information about multiple suppression groups.**

This endpoint will return information for each group ID that you include in your request. To add a group ID to your request, simply append `&id=` followed by the group ID.

Suppressions are a list of email addresses that will not receive content sent under a given [group](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html).

Suppression groups, or [unsubscribe groups](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html), allow you to label a category of content that you regularly send. This gives your recipients the ability to opt out of a specific set of your email. For example, you might define a group for your transactional email, and one for your marketing email so that your users can continue recieving your transactional email witout having to receive your marketing content.

### GET /asm/groups


```php
$query_params = json_decode('{"id": 1}');
$response = $sg->client->asm()->groups()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a suppression group.

**This endpoint allows you to update or change a suppression group.**

Suppression groups, or unsubscribe groups, are specific types or categories of email that you would like your recipients to be able to unsubscribe from. For example: Daily Newsletters, Invoices, System Alerts.

The **name** and **description** of the unsubscribe group will be visible by recipients when they are managing their subscriptions.

Each user can create up to 25 different suppression groups.

### PATCH /asm/groups/{group_id}


```php
$request_body = json_decode('{
  "description": "Suggestions for items our users might like.",
  "id": 103,
  "name": "Item Suggestions"
}');
$group_id = "test_url_param";
$response = $sg->client->asm()->groups()->_($group_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Get information on a single suppression group.

**This endpoint allows you to retrieve a single suppression group.**

Suppression groups, or unsubscribe groups, are specific types or categories of email that you would like your recipients to be able to unsubscribe from. For example: Daily Newsletters, Invoices, System Alerts.

The **name** and **description** of the unsubscribe group will be visible by recipients when they are managing their subscriptions.

Each user can create up to 25 different suppression groups.

### GET /asm/groups/{group_id}


```php
$group_id = "test_url_param";
$response = $sg->client->asm()->groups()->_($group_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a suppression group.

**This endpoint allows you to delete a suppression group.**

You can only delete groups that have not been attached to sent mail in the last 60 days. If a recipient uses the "one-click unsubscribe" option on an email associated with a deleted group, that recipient will be added to the global suppression list.

Suppression groups, or unsubscribe groups, are specific types or categories of email that you would like your recipients to be able to unsubscribe from. For example: Daily Newsletters, Invoices, System Alerts.

The **name** and **description** of the unsubscribe group will be visible by recipients when they are managing their subscriptions.

Each user can create up to 25 different suppression groups.

### DELETE /asm/groups/{group_id}


```php
$group_id = "test_url_param";
$response = $sg->client->asm()->groups()->_($group_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add suppressions to a suppression group

**This endpoint allows you to add email addresses to an unsubscribe group.**

If you attempt to add suppressions to a group that has been deleted or does not exist, the suppressions will be added to the global suppressions list.

Suppressions are recipient email addresses that are added to [unsubscribe groups](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html). Once a recipient's address is on the suppressions list for an unsubscribe group, they will not receive any emails that are tagged with that unsubscribe group.

### POST /asm/groups/{group_id}/suppressions


```php
$request_body = json_decode('{
  "recipient_emails": [
    "test1@example.com",
    "test2@example.com"
  ]
}');
$group_id = "test_url_param";
$response = $sg->client->asm()->groups()->_($group_id)->suppressions()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all suppressions for a suppression group

**This endpoint allows you to retrieve all suppressed email addresses belonging to the given group.**

Suppressions are recipient email addresses that are added to [unsubscribe groups](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html). Once a recipient's address is on the suppressions list for an unsubscribe group, they will not receive any emails that are tagged with that unsubscribe group.

### GET /asm/groups/{group_id}/suppressions


```php
$group_id = "test_url_param";
$response = $sg->client->asm()->groups()->_($group_id)->suppressions()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Search for suppressions within a group

**This endpoint allows you to search a suppression group for multiple suppressions.**

When given a list of email addresses and a group ID, this endpoint will return only the email addresses that have been unsubscribed from the given group.

Suppressions are a list of email addresses that will not receive content sent under a given [group](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html).

### POST /asm/groups/{group_id}/suppressions/search


```php
$request_body = json_decode('{
  "recipient_emails": [
    "exists1@example.com",
    "exists2@example.com",
    "doesnotexists@example.com"
  ]
}');
$group_id = "test_url_param";
$response = $sg->client->asm()->groups()->_($group_id)->suppressions()->search()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a suppression from a suppression group

**This endpoint allows you to remove a suppressed email address from the given suppression group.**

Suppressions are recipient email addresses that are added to [unsubscribe groups](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html). Once a recipient's address is on the suppressions list for an unsubscribe group, they will not receive any emails that are tagged with that unsubscribe group.

### DELETE /asm/groups/{group_id}/suppressions/{email}


```php
$group_id = "test_url_param";
$email = "test_url_param";
$response = $sg->client->asm()->groups()->_($group_id)->suppressions()->_($email)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all suppressions

**This endpoint allows you to retrieve a list of all suppressions.**

Suppressions are a list of email addresses that will not receive content sent under a given [group](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html).

### GET /asm/suppressions


```php
$response = $sg->client->asm()->suppressions()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add recipient addresses to the global suppression group.

**This endpoint allows you to add one or more email addresses to the global suppressions group.**

A global suppression (or global unsubscribe) is an email address of a recipient who does not want to receive any of your messages. A globally suppressed recipient will be removed from any email you send. For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/global_unsubscribes.html).

### POST /asm/suppressions/global


```php
$request_body = json_decode('{
  "recipient_emails": [
    "test1@example.com",
    "test2@example.com"
  ]
}');
$response = $sg->client->asm()->suppressions()->global()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a Global Suppression

**This endpoint allows you to retrieve a global suppression. You can also use this endpoint to confirm if an email address is already globally suppresed.**

If the email address you include in the URL path parameter `{email}` is alreayd globally suppressed, the response will include that email address. If the address you enter for `{email}` is not globally suppressed, an empty JSON object `{}` will be returned.

A global suppression (or global unsubscribe) is an email address of a recipient who does not want to receive any of your messages. A globally suppressed recipient will be removed from any email you send. For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/global_unsubscribes.html).

### GET /asm/suppressions/global/{email}


```php
$email = "test_url_param";
$response = $sg->client->asm()->suppressions()->global()->_($email)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a Global Suppression

**This endpoint allows you to remove an email address from the global suppressions group.**

A global suppression (or global unsubscribe) is an email address of a recipient who does not want to receive any of your messages. A globally suppressed recipient will be removed from any email you send. For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/global_unsubscribes.html).

### DELETE /asm/suppressions/global/{email}


```php
$email = "test_url_param";
$response = $sg->client->asm()->suppressions()->global()->_($email)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all suppression groups for an email address

**This endpoint returns the list of all groups that the given email address has been unsubscribed from.**

Suppressions are a list of email addresses that will not receive content sent under a given [group](https://sendgrid.com/docs/API_Reference/Web_API_v3/Suppression_Management/groups.html).

### GET /asm/suppressions/{email}


```php
$email = "test_url_param";
$response = $sg->client->asm()->suppressions()->_($email)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="browsers"></a>
# BROWSERS

## Retrieve email statistics by browser.

**This endpoint allows you to retrieve your email statistics segmented by browser type.**

**We only store up to 7 days of email activity in our database.** By default, 500 items will be returned per request via the Advanced Stats API endpoints.

Advanced Stats provide a more in-depth view of your email statistics and the actions taken by your recipients. You can segment these statistics by geographic location, device type, client type, browser, and mailbox provider. For more information about statistics, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/index.html).

### GET /browsers/stats


```php
$query_params = json_decode('{"end_date": "2016-04-01", "aggregated_by": "day", "browsers": "test_string", "limit": "test_string", "offset": "test_string", "start_date": "2016-01-01"}');
$response = $sg->client->browsers()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="campaigns"></a>
# CAMPAIGNS

## Create a Campaign

**This endpoint allows you to create a campaign.**

Our Marketing Campaigns API lets you create, manage, send, and schedule campaigns.

Note: In order to send or schedule the campaign, you will be required to provide a subject, sender ID, content (we suggest both html and plain text), and at least one list or segment ID. This information is not required when you create a campaign.

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### POST /campaigns


```php
$request_body = json_decode('{
  "categories": [
    "spring line"
  ],
  "custom_unsubscribe_url": "",
  "html_content": "<html><head><title></title></head><body><p>Check out our spring line!</p></body></html>",
  "ip_pool": "marketing",
  "list_ids": [
    110,
    124
  ],
  "plain_content": "Check out our spring line!",
  "segment_ids": [
    110
  ],
  "sender_id": 124451,
  "subject": "New Products for Spring!",
  "suppression_group_id": 42,
  "title": "March Newsletter"
}');
$response = $sg->client->campaigns()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all Campaigns

**This endpoint allows you to retrieve a list of all of your campaigns.**

Returns campaigns in reverse order they were created (newest first).

Returns an empty array if no campaigns exist.

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### GET /campaigns


```php
$query_params = json_decode('{"limit": 1, "offset": 1}');
$response = $sg->client->campaigns()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a Campaign

Update a campaign. This is especially useful if you only set up the campaign using POST /campaigns, but didn't set many of the parameters.

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### PATCH /campaigns/{campaign_id}


```php
$request_body = json_decode('{
  "categories": [
    "summer line"
  ],
  "html_content": "<html><head><title></title></head><body><p>Check out our summer line!</p></body></html>",
  "plain_content": "Check out our summer line!",
  "subject": "New Products for Summer!",
  "title": "May Newsletter"
}');
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a single campaign

**This endpoint allows you to retrieve a specific campaign.**

Our Marketing Campaigns API lets you create, manage, send, and schedule campaigns.

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### GET /campaigns/{campaign_id}


```php
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a Campaign

**This endpoint allows you to delete a specific campaign.**

Our Marketing Campaigns API lets you create, manage, send, and schedule campaigns.

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### DELETE /campaigns/{campaign_id}


```php
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a Scheduled Campaign

**This endpoint allows to you change the scheduled time and date for a campaign to be sent.**

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### PATCH /campaigns/{campaign_id}/schedules


```php
$request_body = json_decode('{
  "send_at": 1489451436
}');
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->schedules()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Schedule a Campaign

**This endpoint allows you to schedule a specific date and time for your campaign to be sent.**

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### POST /campaigns/{campaign_id}/schedules


```php
$request_body = json_decode('{
  "send_at": 1489771528
}');
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->schedules()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## View Scheduled Time of a Campaign

**This endpoint allows you to retrieve the date and time that the given campaign has been scheduled to be sent.**

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### GET /campaigns/{campaign_id}/schedules


```php
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->schedules()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Unschedule a Scheduled Campaign

**This endpoint allows you to unschedule a campaign that has already been scheduled to be sent.**

A successful unschedule will return a 204.
If the specified campaign is in the process of being sent, the only option is to cancel (a different method).

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### DELETE /campaigns/{campaign_id}/schedules


```php
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->schedules()->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Send a Campaign

**This endpoint allows you to immediately send a campaign at the time you make the API call.**

Normally a POST would have a request body, but since this endpoint is telling us to send a resource that is already created, a request body is not needed.

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### POST /campaigns/{campaign_id}/schedules/now


```php
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->schedules()->now()->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Send a Test Campaign

**This endpoint allows you to send a test campaign.**

To send to multiple addresses, use an array for the JSON "to" value ["one@address","two@address"]

For more information:

* [User Guide > Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html)

### POST /campaigns/{campaign_id}/schedules/test


```php
$request_body = json_decode('{
  "to": "your.email@example.com"
}');
$campaign_id = "test_url_param";
$response = $sg->client->campaigns()->_($campaign_id)->schedules()->test()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="categories"></a>
# CATEGORIES

## Retrieve all categories

**This endpoint allows you to retrieve a list of all of your categories.**

Categories can help organize your email analytics by enabling you to tag emails by type or broad topic. You can define your own custom categories. For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/categories.html).

### GET /categories


```php
$query_params = json_decode('{"category": "test_string", "limit": 1, "offset": 1}');
$response = $sg->client->categories()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve Email Statistics for Categories

**This endpoint allows you to retrieve all of your email statistics for each of your categories.**

If you do not define any query parameters, this endpoint will return a sum for each category in groups of 10.

Categories allow you to group your emails together according to broad topics that you define. For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/categories.html).

### GET /categories/stats


```php
$query_params = json_decode('{"end_date": "2016-04-01", "aggregated_by": "day", "limit": 1, "offset": 1, "start_date": "2016-01-01", "categories": "test_string"}');
$response = $sg->client->categories()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve sums of email stats for each category [Needs: Stats object defined, has category ID?]

**This endpoint allows you to retrieve the total sum of each email statistic for every category over the given date range.**

If you do not define any query parameters, this endpoint will return a sum for each category in groups of 10.

Categories allow you to group your emails together according to broad topics that you define. For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/categories.html).

### GET /categories/stats/sums


```php
$query_params = json_decode('{"end_date": "2016-04-01", "aggregated_by": "day", "limit": 1, "sort_by_metric": "test_string", "offset": 1, "start_date": "2016-01-01", "sort_by_direction": "asc"}');
$response = $sg->client->categories()->stats()->sums()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="clients"></a>
# CLIENTS

## Retrieve email statistics by client type.

**This endpoint allows you to retrieve your email statistics segmented by client type.**

**We only store up to 7 days of email activity in our database.** By default, 500 items will be returned per request via the Advanced Stats API endpoints.

Advanced Stats provide a more in-depth view of your email statistics and the actions taken by your recipients. You can segment these statistics by geographic location, device type, client type, browser, and mailbox provider. For more information about statistics, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/index.html).

### GET /clients/stats


```php
$query_params = json_decode('{"aggregated_by": "day", "start_date": "2016-01-01", "end_date": "2016-04-01"}');
$response = $sg->client->clients()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve stats by a specific client type.

**This endpoint allows you to retrieve your email statistics segmented by a specific client type.**

**We only store up to 7 days of email activity in our database.** By default, 500 items will be returned per request via the Advanced Stats API endpoints.

## Available Client Types
- phone
- tablet
- webmail
- desktop

Advanced Stats provide a more in-depth view of your email statistics and the actions taken by your recipients. You can segment these statistics by geographic location, device type, client type, browser, and mailbox provider. For more information about statistics, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/index.html).

### GET /clients/{client_type}/stats


```php
$query_params = json_decode('{"aggregated_by": "day", "start_date": "2016-01-01", "end_date": "2016-04-01"}');
$client_type = "test_url_param";
$response = $sg->client->clients()->_($client_type)->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="contactdb"></a>
# CONTACTDB

## Create a Custom Field

**This endpoint allows you to create a custom field.**

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### POST /contactdb/custom_fields


```php
$request_body = json_decode('{
  "name": "pet",
  "type": "text"
}');
$response = $sg->client->contactdb()->custom_fields()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all custom fields

**This endpoint allows you to retrieve all custom fields.**

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### GET /contactdb/custom_fields


```php
$response = $sg->client->contactdb()->custom_fields()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a Custom Field

**This endpoint allows you to retrieve a custom field by ID.**

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### GET /contactdb/custom_fields/{custom_field_id}


```php
$custom_field_id = "test_url_param";
$response = $sg->client->contactdb()->custom_fields()->_($custom_field_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a Custom Field

**This endpoint allows you to delete a custom field by ID.**

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### DELETE /contactdb/custom_fields/{custom_field_id}


```php
$custom_field_id = "test_url_param";
$response = $sg->client->contactdb()->custom_fields()->_($custom_field_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create a List

**This endpoint allows you to create a list for your recipients.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### POST /contactdb/lists


```php
$request_body = json_decode('{
  "name": "your list name"
}');
$response = $sg->client->contactdb()->lists()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all lists

**This endpoint allows you to retrieve all of your recipient lists. If you don't have any lists, an empty array will be returned.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### GET /contactdb/lists


```php
$response = $sg->client->contactdb()->lists()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete Multiple lists

**This endpoint allows you to delete multiple recipient lists.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### DELETE /contactdb/lists


```php
$request_body = json_decode('[
  1,
  2,
  3,
  4
]');
$response = $sg->client->contactdb()->lists()->delete($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a List

**This endpoint allows you to update the name of one of your recipient lists.**


The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### PATCH /contactdb/lists/{list_id}


```php
$request_body = json_decode('{
  "name": "newlistname"
}');
$query_params = json_decode('{"list_id": 1}');
$list_id = "test_url_param";
$response = $sg->client->contactdb()->lists()->_($list_id)->patch($request_body, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a single list

This endpoint allows you to retrieve a single recipient list.

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### GET /contactdb/lists/{list_id}


```php
$query_params = json_decode('{"list_id": 1}');
$list_id = "test_url_param";
$response = $sg->client->contactdb()->lists()->_($list_id)->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a List

**This endpoint allows you to delete a specific recipient list with the given ID.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### DELETE /contactdb/lists/{list_id}


```php
$query_params = json_decode('{"delete_contacts": "true"}');
$list_id = "test_url_param";
$response = $sg->client->contactdb()->lists()->_($list_id)->delete(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add Multiple Recipients to a List

**This endpoint allows you to add multiple recipients to a list.**

Adds existing recipients to a list, passing in the recipient IDs to add. Recipient IDs should be passed exactly as they are returned from recipient endpoints.

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### POST /contactdb/lists/{list_id}/recipients


```php
$request_body = json_decode('[
  "recipient_id1",
  "recipient_id2"
]');
$list_id = "test_url_param";
$response = $sg->client->contactdb()->lists()->_($list_id)->recipients()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all recipients on a List

**This endpoint allows you to retrieve all recipients on the list with the given ID.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### GET /contactdb/lists/{list_id}/recipients


```php
$query_params = json_decode('{"page": 1, "page_size": 1, "list_id": 1}');
$list_id = "test_url_param";
$response = $sg->client->contactdb()->lists()->_($list_id)->recipients()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add a Single Recipient to a List

**This endpoint allows you to add a single recipient to a list.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### POST /contactdb/lists/{list_id}/recipients/{recipient_id}


```php
$list_id = "test_url_param";
$recipient_id = "test_url_param";
$response = $sg->client->contactdb()->lists()->_($list_id)->recipients()->_($recipient_id)->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a Single Recipient from a Single List

**This endpoint allows you to delete a single recipient from a list.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### DELETE /contactdb/lists/{list_id}/recipients/{recipient_id}


```php
$query_params = json_decode('{"recipient_id": 1, "list_id": 1}');
$list_id = "test_url_param";
$recipient_id = "test_url_param";
$response = $sg->client->contactdb()->lists()->_($list_id)->recipients()->_($recipient_id)->delete(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Recipient

**This endpoint allows you to update one or more recipients.**

The body of an API call to this endpoint must include an array of one or more recipient objects.

It is of note that you can add custom field data as parameters on recipient objects. We have provided an example using some of the default custom fields SendGrid provides.

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### PATCH /contactdb/recipients


```php
$request_body = json_decode('[
  {
    "email": "jones@example.com",
    "first_name": "Guy",
    "last_name": "Jones"
  }
]');
$response = $sg->client->contactdb()->recipients()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add recipients

**This endpoint allows you to add a Marketing Campaigns recipient.**

It is of note that you can add custom field data as a parameter on this endpoint. We have provided an example using some of the default custom fields SendGrid provides.

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### POST /contactdb/recipients


```php
$request_body = json_decode('[
  {
    "age": 25,
    "email": "example@example.com",
    "first_name": "",
    "last_name": "User"
  },
  {
    "age": 25,
    "email": "example2@example.com",
    "first_name": "Example",
    "last_name": "User"
  }
]');
$response = $sg->client->contactdb()->recipients()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve recipients

**This endpoint allows you to retrieve all of your Marketing Campaigns recipients.**

Batch deletion of a page makes it possible to receive an empty page of recipients before reaching the end of
the list of recipients. To avoid this issue; iterate over pages until a 404 is retrieved.

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### GET /contactdb/recipients


```php
$query_params = json_decode('{"page": 1, "page_size": 1}');
$response = $sg->client->contactdb()->recipients()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete Recipient

**This endpoint allows you to deletes one or more recipients.**

The body of an API call to this endpoint must include an array of recipient IDs of the recipients you want to delete.

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### DELETE /contactdb/recipients


```php
$request_body = json_decode('[
  "recipient_id1",
  "recipient_id2"
]');
$response = $sg->client->contactdb()->recipients()->delete($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve the count of billable recipients

**This endpoint allows you to retrieve the number of Marketing Campaigns recipients that you will be billed for.**

You are billed for marketing campaigns based on the highest number of recipients you have had in your account at one time. This endpoint will allow you to know the current billable count value.

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### GET /contactdb/recipients/billable_count


```php
$response = $sg->client->contactdb()->recipients()->billable_count()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a Count of Recipients

**This endpoint allows you to retrieve the total number of Marketing Campaigns recipients.**

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### GET /contactdb/recipients/count


```php
$response = $sg->client->contactdb()->recipients()->count()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve recipients matching search criteria

**This endpoint allows you to perform a search on all of your Marketing Campaigns recipients.**

field_name:

* is a variable that is substituted for your actual custom field name from your recipient.
* Text fields must be url-encoded. Date fields are searchable only by unix timestamp (e.g. 2/2/2015 becomes 1422835200)
* If field_name is a 'reserved' date field, such as created_at or updated_at, the system will internally convert
your epoch time to a date range encompassing the entire day. For example, an epoch time of 1422835600 converts to
Mon, 02 Feb 2015 00:06:40 GMT, but internally the system will search from Mon, 02 Feb 2015 00:00:00 GMT through
Mon, 02 Feb 2015 23:59:59 GMT.

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### GET /contactdb/recipients/search


```php
$query_params = json_decode('{"{field_name}": "test_string"}');
$response = $sg->client->contactdb()->recipients()->search()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a single recipient

**This endpoint allows you to retrieve a single recipient by ID from your contact database.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### GET /contactdb/recipients/{recipient_id}


```php
$recipient_id = "test_url_param";
$response = $sg->client->contactdb()->recipients()->_($recipient_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a Recipient

**This endpoint allows you to delete a single recipient with the given ID from your contact database.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### DELETE /contactdb/recipients/{recipient_id}


```php
$recipient_id = "test_url_param";
$response = $sg->client->contactdb()->recipients()->_($recipient_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve the lists that a recipient is on

**This endpoint allows you to retrieve the lists that a given recipient belongs to.**

Each recipient can be on many lists. This endpoint gives you all of the lists that any one recipient has been added to.

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

### GET /contactdb/recipients/{recipient_id}/lists


```php
$recipient_id = "test_url_param";
$response = $sg->client->contactdb()->recipients()->_($recipient_id)->lists()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve reserved fields

**This endpoint allows you to list all fields that are reserved and can't be used for custom field names.**

The contactdb is a database of your contacts for [SendGrid Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html).

### GET /contactdb/reserved_fields


```php
$response = $sg->client->contactdb()->reserved_fields()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create a Segment

**This endpoint allows you to create a segment.**

All recipients in your contactdb will be added or removed automatically depending on whether they match the criteria for this segment.

List Id:

* Send this to segment from an existing list
* Don't send this in order to segment from your entire contactdb.

Valid operators for create and update depend on the type of the field you are segmenting:

* **Dates:** "eq", "ne", "lt" (before), "gt" (after)
* **Text:** "contains", "eq" (is - matches the full field), "ne" (is not - matches any field where the entire field is not the condition value)
* **Numbers:** "eq", "lt", "gt"
* **Email Clicks and Opens:** "eq" (opened), "ne" (not opened)

Segment conditions using "eq" or "ne" for email clicks and opens should provide a "field" of either *clicks.campaign_identifier* or *opens.campaign_identifier*. The condition value should be a string containing the id of a completed campaign.

Segments may contain multiple condtions, joined by an "and" or "or" in the "and_or" field. The first condition in the conditions list must have an empty "and_or", and subsequent conditions must all specify an "and_or".

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

For more information about segments in Marketing Campaigns, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/lists.html#-Create-a-Segment).

### POST /contactdb/segments


```php
$request_body = json_decode('{
  "conditions": [
    {
      "and_or": "",
      "field": "last_name",
      "operator": "eq",
      "value": "Miller"
    },
    {
      "and_or": "and",
      "field": "last_clicked",
      "operator": "gt",
      "value": "01/02/2015"
    },
    {
      "and_or": "or",
      "field": "clicks.campaign_identifier",
      "operator": "eq",
      "value": "513"
    }
  ],
  "list_id": 4,
  "name": "Last Name Miller"
}');
$response = $sg->client->contactdb()->segments()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all segments

**This endpoint allows you to retrieve all of your segments.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

For more information about segments in Marketing Campaigns, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/lists.html#-Create-a-Segment).

### GET /contactdb/segments


```php
$response = $sg->client->contactdb()->segments()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a segment

**This endpoint allows you to update a segment.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

For more information about segments in Marketing Campaigns, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/lists.html#-Create-a-Segment).

### PATCH /contactdb/segments/{segment_id}


```php
$request_body = json_decode('{
  "conditions": [
    {
      "and_or": "",
      "field": "last_name",
      "operator": "eq",
      "value": "Miller"
    }
  ],
  "list_id": 5,
  "name": "The Millers"
}');
$query_params = json_decode('{"segment_id": "test_string"}');
$segment_id = "test_url_param";
$response = $sg->client->contactdb()->segments()->_($segment_id)->patch($request_body, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a segment

**This endpoint allows you to retrieve a single segment with the given ID.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

For more information about segments in Marketing Campaigns, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/lists.html#-Create-a-Segment).

### GET /contactdb/segments/{segment_id}


```php
$query_params = json_decode('{"segment_id": 1}');
$segment_id = "test_url_param";
$response = $sg->client->contactdb()->segments()->_($segment_id)->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a segment

**This endpoint allows you to delete a segment from your recipients database.**

You also have the option to delete all the contacts from your Marketing Campaigns recipient database who were in this segment.

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

For more information about segments in Marketing Campaigns, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/lists.html#-Create-a-Segment).

### DELETE /contactdb/segments/{segment_id}


```php
$query_params = json_decode('{"delete_contacts": "true"}');
$segment_id = "test_url_param";
$response = $sg->client->contactdb()->segments()->_($segment_id)->delete(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve recipients on a segment

**This endpoint allows you to retrieve all of the recipients in a segment with the given ID.**

The Contacts API helps you manage your [Marketing Campaigns](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/index.html) recipients.

For more information about segments in Marketing Campaigns, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/lists.html#-Create-a-Segment).

### GET /contactdb/segments/{segment_id}/recipients


```php
$query_params = json_decode('{"page": 1, "page_size": 1}');
$segment_id = "test_url_param";
$response = $sg->client->contactdb()->segments()->_($segment_id)->recipients()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="devices"></a>
# DEVICES

## Retrieve email statistics by device type.

**This endpoint allows you to retrieve your email statistics segmented by the device type.**

**We only store up to 7 days of email activity in our database.** By default, 500 items will be returned per request via the Advanced Stats API endpoints.

## Available Device Types
| **Device** | **Description** | **Example** |
|---|---|---|
| Desktop | Email software on desktop computer. | I.E., Outlook, Sparrow, or Apple Mail. |
| Webmail |	A web-based email client. | I.E., Yahoo, Google, AOL, or Outlook.com. |
| Phone | A smart phone. | iPhone, Android, Blackberry, etc.
| Tablet | A tablet computer. | iPad, android based tablet, etc. |
| Other | An unrecognized device. |

Advanced Stats provide a more in-depth view of your email statistics and the actions taken by your recipients. You can segment these statistics by geographic location, device type, client type, browser, and mailbox provider. For more information about statistics, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/index.html).

### GET /devices/stats


```php
$query_params = json_decode('{"aggregated_by": "day", "limit": 1, "start_date": "2016-01-01", "end_date": "2016-04-01", "offset": 1}');
$response = $sg->client->devices()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="geo"></a>
# GEO

## Retrieve email statistics by country and state/province.

**This endpoint allows you to retrieve your email statistics segmented by country and state/province.**

**We only store up to 7 days of email activity in our database.** By default, 500 items will be returned per request via the Advanced Stats API endpoints.

Advanced Stats provide a more in-depth view of your email statistics and the actions taken by your recipients. You can segment these statistics by geographic location, device type, client type, browser, and mailbox provider. For more information about statistics, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/index.html).

### GET /geo/stats


```php
$query_params = json_decode('{"end_date": "2016-04-01", "country": "US", "aggregated_by": "day", "limit": 1, "offset": 1, "start_date": "2016-01-01"}');
$response = $sg->client->geo()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="ips"></a>
# IPS

## Retrieve all IP addresses

**This endpoint allows you to retrieve a list of all assigned and unassigned IPs.**

Response includes warm up status, pools, assigned subusers, and whitelabel info. The start_date field corresponds to when warmup started for that IP.

A single IP address or a range of IP addresses may be dedicated to an account in order to send email for multiple domains. The reputation of this IP is based on the aggregate performance of all the senders who use it.

### GET /ips


```php
$query_params = json_decode('{"subuser": "test_string", "ip": "test_string", "limit": 1, "exclude_whitelabels": "true", "offset": 1}');
$response = $sg->client->ips()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all assigned IPs

**This endpoint allows you to retrieve only assigned IP addresses.**

A single IP address or a range of IP addresses may be dedicated to an account in order to send email for multiple domains. The reputation of this IP is based on the aggregate performance of all the senders who use it.

### GET /ips/assigned


```php
$response = $sg->client->ips()->assigned()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create an IP pool.

**This endpoint allows you to create an IP pool.**

**Each user can create up to 10 different IP pools.**

IP Pools allow you to group your dedicated SendGrid IP addresses together. For example, you could create separate pools for your transactional and marketing email. When sending marketing emails, specify that you want to use the marketing IP pool. This allows you to maintain separate reputations for your different email traffic.

IP pools can only be used with whitelabeled IP addresses.

If an IP pool is NOT specified for an email, it will use any IP available, including ones in pools.

### POST /ips/pools


```php
$request_body = json_decode('{
  "name": "marketing"
}');
$response = $sg->client->ips()->pools()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all IP pools.

**This endpoint allows you to retreive all of your IP pools.**

IP Pools allow you to group your dedicated SendGrid IP addresses together. For example, you could create separate pools for your transactional and marketing email. When sending marketing emails, specify that you want to use the marketing IP pool. This allows you to maintain separate reputations for your different email traffic.

IP pools can only be used with whitelabeled IP addresses.

If an IP pool is NOT specified for an email, it will use any IP available, including ones in pools.

### GET /ips/pools


```php
$response = $sg->client->ips()->pools()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update an IP pools name.

**This endpoint allows you to update the name of an IP pool.**

IP Pools allow you to group your dedicated SendGrid IP addresses together. For example, you could create separate pools for your transactional and marketing email. When sending marketing emails, specify that you want to use the marketing IP pool. This allows you to maintain separate reputations for your different email traffic.

IP pools can only be used with whitelabeled IP addresses.

If an IP pool is NOT specified for an email, it will use any IP available, including ones in pools.

### PUT /ips/pools/{pool_name}


```php
$request_body = json_decode('{
  "name": "new_pool_name"
}');
$pool_name = "test_url_param";
$response = $sg->client->ips()->pools()->_($pool_name)->put($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all IPs in a specified pool.

**This endpoint allows you to list all of the IP addresses that are in a specific IP pool.**

IP Pools allow you to group your dedicated SendGrid IP addresses together. For example, you could create separate pools for your transactional and marketing email. When sending marketing emails, specify that you want to use the marketing IP pool. This allows you to maintain separate reputations for your different email traffic.

IP pools can only be used with whitelabeled IP addresses.

If an IP pool is NOT specified for an email, it will use any IP available, including ones in pools.

### GET /ips/pools/{pool_name}


```php
$pool_name = "test_url_param";
$response = $sg->client->ips()->pools()->_($pool_name)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete an IP pool.

**This endpoint allows you to delete an IP pool.**

IP Pools allow you to group your dedicated SendGrid IP addresses together. For example, you could create separate pools for your transactional and marketing email. When sending marketing emails, specify that you want to use the marketing IP pool. This allows you to maintain separate reputations for your different email traffic.

IP pools can only be used with whitelabeled IP addresses.

If an IP pool is NOT specified for an email, it will use any IP available, including ones in pools.

### DELETE /ips/pools/{pool_name}


```php
$pool_name = "test_url_param";
$response = $sg->client->ips()->pools()->_($pool_name)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add an IP address to a pool

**This endpoint allows you to add an IP address to an IP pool.**

You can add the same IP address to multiple pools. It may take up to 60 seconds for your IP address to be added to a pool after your request is made.

A single IP address or a range of IP addresses may be dedicated to an account in order to send email for multiple domains. The reputation of this IP is based on the aggregate performance of all the senders who use it.

### POST /ips/pools/{pool_name}/ips


```php
$request_body = json_decode('{
  "ip": "0.0.0.0"
}');
$pool_name = "test_url_param";
$response = $sg->client->ips()->pools()->_($pool_name)->ips()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Remove an IP address from a pool.

**This endpoint allows you to remove an IP address from an IP pool.**

The same IP address can be added to multiple IP pools.

A single IP address or a range of IP addresses may be dedicated to an account in order to send email for multiple domains. The reputation of this IP is based on the aggregate performance of all the senders who use it.

### DELETE /ips/pools/{pool_name}/ips/{ip}


```php
$pool_name = "test_url_param";
$ip = "test_url_param";
$response = $sg->client->ips()->pools()->_($pool_name)->ips()->_($ip)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add an IP to warmup

**This endpoint allows you to enter an IP address into warmup mode.**

SendGrid can automatically warm up dedicated IP addresses by limiting the amount of mail that can be sent through them per hour, with the limit determined by how long the IP address has been in warmup. See the [warmup schedule](https://sendgrid.com/docs/API_Reference/Web_API_v3/IP_Management/ip_warmup_schedule.html) for more details on how SendGrid limits your email traffic for IPs in warmup.

For more general information about warming up IPs, please see our [Classroom](https://sendgrid.com/docs/Classroom/Deliver/Delivery_Introduction/warming_up_ips.html).

### POST /ips/warmup


```php
$request_body = json_decode('{
  "ip": "0.0.0.0"
}');
$response = $sg->client->ips()->warmup()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all IPs currently in warmup

**This endpoint allows you to retrieve all of your IP addresses that are currently warming up.**

SendGrid can automatically warm up dedicated IP addresses by limiting the amount of mail that can be sent through them per hour, with the limit determined by how long the IP address has been in warmup. See the [warmup schedule](https://sendgrid.com/docs/API_Reference/Web_API_v3/IP_Management/ip_warmup_schedule.html) for more details on how SendGrid limits your email traffic for IPs in warmup.

For more general information about warming up IPs, please see our [Classroom](https://sendgrid.com/docs/Classroom/Deliver/Delivery_Introduction/warming_up_ips.html).

### GET /ips/warmup


```php
$response = $sg->client->ips()->warmup()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve warmup status for a specific IP address

**This endpoint allows you to retrieve the warmup status for a specific IP address.**

SendGrid can automatically warm up dedicated IP addresses by limiting the amount of mail that can be sent through them per hour, with the limit determined by how long the IP address has been in warmup. See the [warmup schedule](https://sendgrid.com/docs/API_Reference/Web_API_v3/IP_Management/ip_warmup_schedule.html) for more details on how SendGrid limits your email traffic for IPs in warmup.

For more general information about warming up IPs, please see our [Classroom](https://sendgrid.com/docs/Classroom/Deliver/Delivery_Introduction/warming_up_ips.html).

### GET /ips/warmup/{ip_address}


```php
$ip_address = "test_url_param";
$response = $sg->client->ips()->warmup()->_($ip_address)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Remove an IP from warmup

**This endpoint allows you to remove an IP address from warmup mode.**

SendGrid can automatically warm up dedicated IP addresses by limiting the amount of mail that can be sent through them per hour, with the limit determined by how long the IP address has been in warmup. See the [warmup schedule](https://sendgrid.com/docs/API_Reference/Web_API_v3/IP_Management/ip_warmup_schedule.html) for more details on how SendGrid limits your email traffic for IPs in warmup.

For more general information about warming up IPs, please see our [Classroom](https://sendgrid.com/docs/Classroom/Deliver/Delivery_Introduction/warming_up_ips.html).

### DELETE /ips/warmup/{ip_address}


```php
$ip_address = "test_url_param";
$response = $sg->client->ips()->warmup()->_($ip_address)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all IP pools an IP address belongs to

**This endpoint allows you to see which IP pools a particular IP address has been added to.**

The same IP address can be added to multiple IP pools.

A single IP address or a range of IP addresses may be dedicated to an account in order to send email for multiple domains. The reputation of this IP is based on the aggregate performance of all the senders who use it.

### GET /ips/{ip_address}


```php
$ip_address = "test_url_param";
$response = $sg->client->ips()->_($ip_address)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="mail"></a>
# MAIL

## Create a batch ID

**This endpoint allows you to generate a new batch ID. This batch ID can be associated with scheduled sends via the mail/send endpoint.**

If you set the SMTPAPI header `batch_id`, it allows you to then associate multiple scheduled mail/send requests together with the same ID. Then at anytime up to 10 minutes before the schedule date, you can cancel all of the mail/send requests that have this batch ID by calling the Cancel Scheduled Send endpoint.

More Information:

* [Scheduling Parameters > Batch ID](https://sendgrid.com/docs/API_Reference/SMTP_API/scheduling_parameters.html)

### POST /mail/batch


```php
$response = $sg->client->mail()->batch()->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Validate batch ID

**This endpoint allows you to validate a batch ID.**

If you set the SMTPAPI header `batch_id`, it allows you to then associate multiple scheduled mail/send requests together with the same ID. Then at anytime up to 10 minutes before the schedule date, you can cancel all of the mail/send requests that have this batch ID by calling the Cancel Scheduled Send endpoint.

More Information:

* [Scheduling Parameters > Batch ID](https://sendgrid.com/docs/API_Reference/SMTP_API/scheduling_parameters.html)

### GET /mail/batch/{batch_id}


```php
$batch_id = "test_url_param";
$response = $sg->client->mail()->batch()->_($batch_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## v3 Mail Send

This endpoint allows you to send email over SendGrids v3 Web API, the most recent version of our API. If you are looking for documentation about the v2 Mail Send endpoint, please see our [v2 API Reference](https://sendgrid.com/docs/API_Reference/Web_API/mail.html).

* Top level parameters are referred to as "global".
* Individual fields within the personalizations array will override any other global, or message level, parameters that are defined outside of personalizations.

For an overview of the v3 Mail Send endpoint, please visit our [v3 API Reference](https://sendgrid.com/docs/API_Reference/Web_API_v3/Mail/index.html)

For more detailed information about how to use the v3 Mail Send endpoint, please visit our [Classroom](https://sendgrid.com/docs/Classroom/Send/v3_Mail_Send/index.html).

### POST /mail/send

This endpoint has a helper, check it out [here](https://github.com/sendgrid/sendgrid-php/blob/master/lib/helpers/mail/README.md).

```php
$request_body = json_decode('{
  "asm": {
    "group_id": 1,
    "groups_to_display": [
      1,
      2,
      3
    ]
  },
  "attachments": [
    {
      "content": "[BASE64 encoded content block here]",
      "content_id": "ii_139db99fdb5c3704",
      "disposition": "inline",
      "filename": "file1.jpg",
      "name": "file1",
      "type": "jpg"
    }
  ],
  "batch_id": "[YOUR BATCH ID GOES HERE]",
  "categories": [
    "category1",
    "category2"
  ],
  "content": [
    {
      "type": "text/html",
      "value": "<html><p>Hello, world!</p><img src=[CID GOES HERE]></img></html>"
    }
  ],
  "custom_args": {
    "New Argument 1": "New Value 1",
    "activationAttempt": "1",
    "customerAccountNumber": "[CUSTOMER ACCOUNT NUMBER GOES HERE]"
  },
  "from": {
    "email": "sam.smith@example.com",
    "name": "Sam Smith"
  },
  "headers": {},
  "ip_pool_name": "[YOUR POOL NAME GOES HERE]",
  "mail_settings": {
    "bcc": {
      "email": "ben.doe@example.com",
      "enable": true
    },
    "bypass_list_management": {
      "enable": true
    },
    "footer": {
      "enable": true,
      "html": "<p>Thanks</br>The SendGrid Team</p>",
      "text": "Thanks,/n The SendGrid Team"
    },
    "sandbox_mode": {
      "enable": false
    },
    "spam_check": {
      "enable": true,
      "post_to_url": "http://example.com/compliance",
      "threshold": 3
    }
  },
  "personalizations": [
    {
      "bcc": [
        {
          "email": "sam.doe@example.com",
          "name": "Sam Doe"
        }
      ],
      "cc": [
        {
          "email": "jane.doe@example.com",
          "name": "Jane Doe"
        }
      ],
      "custom_args": {
        "New Argument 1": "New Value 1",
        "activationAttempt": "1",
        "customerAccountNumber": "[CUSTOMER ACCOUNT NUMBER GOES HERE]"
      },
      "headers": {
        "X-Accept-Language": "en",
        "X-Mailer": "MyApp"
      },
      "send_at": 1409348513,
      "subject": "Hello, World!",
      "substitutions": {
        "id": "substitutions",
        "type": "object"
      },
      "to": [
        {
          "email": "john.doe@example.com",
          "name": "John Doe"
        }
      ]
    }
  ],
  "reply_to": {
    "email": "sam.smith@example.com",
    "name": "Sam Smith"
  },
  "sections": {
    "section": {
      ":sectionName1": "section 1 text",
      ":sectionName2": "section 2 text"
    }
  },
  "send_at": 1409348513,
  "subject": "Hello, World!",
  "template_id": "[YOUR TEMPLATE ID GOES HERE]",
  "tracking_settings": {
    "click_tracking": {
      "enable": true,
      "enable_text": true
    },
    "ganalytics": {
      "enable": true,
      "utm_campaign": "[NAME OF YOUR REFERRER SOURCE]",
      "utm_content": "[USE THIS SPACE TO DIFFERENTIATE YOUR EMAIL FROM ADS]",
      "utm_medium": "[NAME OF YOUR MARKETING MEDIUM e.g. email]",
      "utm_name": "[NAME OF YOUR CAMPAIGN]",
      "utm_term": "[IDENTIFY PAID KEYWORDS HERE]"
    },
    "open_tracking": {
      "enable": true,
      "substitution_tag": "%opentrack"
    },
    "subscription_tracking": {
      "enable": true,
      "html": "If you would like to unsubscribe and stop receiving these emails <% clickhere %>.",
      "substitution_tag": "<%click here%>",
      "text": "If you would like to unsubscribe and stop receiveing these emails <% click here %>."
    }
  }
}');
$response = $sg->client->mail()->send()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="mail_settings"></a>
# MAIL SETTINGS

## Retrieve all mail settings

**This endpoint allows you to retrieve a list of all mail settings.**

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings


```php
$query_params = json_decode('{"limit": 1, "offset": 1}');
$response = $sg->client->mail_settings()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update address whitelist mail settings

**This endpoint allows you to update your current email address whitelist settings.**

The address whitelist setting whitelists a specified email address or domain for which mail should never be suppressed. For example, you own the domain example.com, and one or more of your recipients use email@example.com addresses, by placing example.com in the address whitelist setting, all bounces, blocks, and unsubscribes logged for that domain will be ignored and sent as if under normal sending conditions.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/address_whitelist


```php
$request_body = json_decode('{
  "enabled": true,
  "list": [
    "email1@example.com",
    "example.com"
  ]
}');
$response = $sg->client->mail_settings()->address_whitelist()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve address whitelist mail settings

**This endpoint allows you to retrieve your current email address whitelist settings.**

The address whitelist setting whitelists a specified email address or domain for which mail should never be suppressed. For example, you own the domain example.com, and one or more of your recipients use email@example.com addresses, by placing example.com in the address whitelist setting, all bounces, blocks, and unsubscribes logged for that domain will be ignored and sent as if under normal sending conditions.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/address_whitelist


```php
$response = $sg->client->mail_settings()->address_whitelist()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update BCC mail settings

**This endpoint allows you to update your current BCC mail settings.**

When the BCC mail setting is enabled, SendGrid will automatically send a blind carbon copy (BCC) to an address for every email sent without adding that address to the header. Please note that only one email address may be entered in this field, if you wish to distribute BCCs to multiple addresses you will need to create a distribution group or use forwarding rules.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/bcc


```php
$request_body = json_decode('{
  "email": "email@example.com",
  "enabled": false
}');
$response = $sg->client->mail_settings()->bcc()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all BCC mail settings

**This endpoint allows you to retrieve your current BCC mail settings.**

When the BCC mail setting is enabled, SendGrid will automatically send a blind carbon copy (BCC) to an address for every email sent without adding that address to the header. Please note that only one email address may be entered in this field, if you wish to distribute BCCs to multiple addresses you will need to create a distribution group or use forwarding rules.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/bcc


```php
$response = $sg->client->mail_settings()->bcc()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update bounce purge mail settings

**This endpoint allows you to update your current bounce purge settings.**

This setting allows you to set a schedule for SendGrid to automatically delete contacts from your soft and hard bounce suppression lists.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/bounce_purge


```php
$request_body = json_decode('{
  "enabled": true,
  "hard_bounces": 5,
  "soft_bounces": 5
}');
$response = $sg->client->mail_settings()->bounce_purge()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve bounce purge mail settings

**This endpoint allows you to retrieve your current bounce purge settings.**

This setting allows you to set a schedule for SendGrid to automatically delete contacts from your soft and hard bounce suppression lists.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/bounce_purge


```php
$response = $sg->client->mail_settings()->bounce_purge()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update footer mail settings

**This endpoint allows you to update your current Footer mail settings.**

The footer setting will insert a custom footer at the bottom of the text and HTML bodies. Use the embedded HTML editor and plain text entry fields to create the content of the footers to be inserted into your emails.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/footer


```php
$request_body = json_decode('{
  "enabled": true,
  "html_content": "...",
  "plain_content": "..."
}');
$response = $sg->client->mail_settings()->footer()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve footer mail settings

**This endpoint allows you to retrieve your current Footer mail settings.**

The footer setting will insert a custom footer at the bottom of the text and HTML bodies. Use the embedded HTML editor and plain text entry fields to create the content of the footers to be inserted into your emails.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/footer


```php
$response = $sg->client->mail_settings()->footer()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update forward bounce mail settings

**This endpoint allows you to update your current bounce forwarding mail settings.**

Activating this setting allows you to specify an email address to which bounce reports are forwarded.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/forward_bounce


```php
$request_body = json_decode('{
  "email": "example@example.com",
  "enabled": true
}');
$response = $sg->client->mail_settings()->forward_bounce()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve forward bounce mail settings

**This endpoint allows you to retrieve your current bounce forwarding mail settings.**

Activating this setting allows you to specify an email address to which bounce reports are forwarded.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/forward_bounce


```php
$response = $sg->client->mail_settings()->forward_bounce()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update forward spam mail settings

**This endpoint allows you to update your current Forward Spam mail settings.**

Enabling the forward spam setting allows you to specify an email address to which spam reports will be forwarded.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/forward_spam


```php
$request_body = json_decode('{
  "email": "",
  "enabled": false
}');
$response = $sg->client->mail_settings()->forward_spam()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve forward spam mail settings

**This endpoint allows you to retrieve your current Forward Spam mail settings.**

Enabling the forward spam setting allows you to specify an email address to which spam reports will be forwarded.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/forward_spam


```php
$response = $sg->client->mail_settings()->forward_spam()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update plain content mail settings

**This endpoint allows you to update your current Plain Content mail settings.**

The plain content setting will automatically convert any plain text emails that you send to HTML before sending.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/plain_content


```php
$request_body = json_decode('{
  "enabled": false
}');
$response = $sg->client->mail_settings()->plain_content()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve plain content mail settings

**This endpoint allows you to retrieve your current Plain Content mail settings.**

The plain content setting will automatically convert any plain text emails that you send to HTML before sending.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/plain_content


```php
$response = $sg->client->mail_settings()->plain_content()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update spam check mail settings

**This endpoint allows you to update your current spam checker mail settings.**

The spam checker filter notifies you when emails are detected that exceed a predefined spam threshold.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/spam_check


```php
$request_body = json_decode('{
  "enabled": true,
  "max_score": 5,
  "url": "url"
}');
$response = $sg->client->mail_settings()->spam_check()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve spam check mail settings

**This endpoint allows you to retrieve your current Spam Checker mail settings.**

The spam checker filter notifies you when emails are detected that exceed a predefined spam threshold.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/spam_check


```php
$response = $sg->client->mail_settings()->spam_check()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update template mail settings

**This endpoint allows you to update your current legacy email template settings.**

This setting refers to our original email templates. We currently support more fully featured [transactional templates](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

The legacy email template setting wraps an HTML template around your email content. This can be useful for sending out marketing email and/or other HTML formatted messages.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### PATCH /mail_settings/template


```php
$request_body = json_decode('{
  "enabled": true,
  "html_content": "<% body %>"
}');
$response = $sg->client->mail_settings()->template()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve legacy template mail settings

**This endpoint allows you to retrieve your current legacy email template settings.**

This setting refers to our original email templates. We currently support more fully featured [transactional templates](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

The legacy email template setting wraps an HTML template around your email content. This can be useful for sending out marketing email and/or other HTML formatted messages.

Mail settings allow you to tell SendGrid specific things to do to every email that you send to your recipients over SendGrids [Web API](https://sendgrid.com/docs/API_Reference/Web_API/mail.html) or [SMTP Relay](https://sendgrid.com/docs/API_Reference/SMTP_API/index.html).

### GET /mail_settings/template


```php
$response = $sg->client->mail_settings()->template()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="mailbox_providers"></a>
# MAILBOX PROVIDERS

## Retrieve email statistics by mailbox provider.

**This endpoint allows you to retrieve your email statistics segmented by recipient mailbox provider.**

**We only store up to 7 days of email activity in our database.** By default, 500 items will be returned per request via the Advanced Stats API endpoints.

Advanced Stats provide a more in-depth view of your email statistics and the actions taken by your recipients. You can segment these statistics by geographic location, device type, client type, browser, and mailbox provider. For more information about statistics, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/index.html).

### GET /mailbox_providers/stats


```php
$query_params = json_decode('{"end_date": "2016-04-01", "mailbox_providers": "test_string", "aggregated_by": "day", "limit": 1, "offset": 1, "start_date": "2016-01-01"}');
$response = $sg->client->mailbox_providers()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="partner_settings"></a>
# PARTNER SETTINGS

## Returns a list of all partner settings.

**This endpoint allows you to retrieve a list of all partner settings that you can enable.**

Our partner settings allow you to integrate your SendGrid account with our partners to increase your SendGrid experience and functionality. For more information about our partners, and how you can begin integrating with them, please visit our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/partners.html).

### GET /partner_settings


```php
$query_params = json_decode('{"limit": 1, "offset": 1}');
$response = $sg->client->partner_settings()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Updates New Relic partner settings.

**This endpoint allows you to update or change your New Relic partner settings.**

Our partner settings allow you to integrate your SendGrid account with our partners to increase your SendGrid experience and functionality. For more information about our partners, and how you can begin integrating with them, please visit our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/partners.html).

By integrating with New Relic, you can send your SendGrid email statistics to your New Relic Dashboard. If you enable this setting, your stats will be sent to New Relic every 5 minutes. You will need your New Relic License Key to enable this setting. For more information, please see our [Classroom](https://sendgrid.com/docs/Classroom/Track/Collecting_Data/new_relic.html).

### PATCH /partner_settings/new_relic


```php
$request_body = json_decode('{
  "enable_subuser_statistics": true,
  "enabled": true,
  "license_key": ""
}');
$response = $sg->client->partner_settings()->new_relic()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Returns all New Relic partner settings.

**This endpoint allows you to retrieve your current New Relic partner settings.**

Our partner settings allow you to integrate your SendGrid account with our partners to increase your SendGrid experience and functionality. For more information about our partners, and how you can begin integrating with them, please visit our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/partners.html).

By integrating with New Relic, you can send your SendGrid email statistics to your New Relic Dashboard. If you enable this setting, your stats will be sent to New Relic every 5 minutes. You will need your New Relic License Key to enable this setting. For more information, please see our [Classroom](https://sendgrid.com/docs/Classroom/Track/Collecting_Data/new_relic.html).

### GET /partner_settings/new_relic


```php
$response = $sg->client->partner_settings()->new_relic()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="scopes"></a>
# SCOPES

## Retrieve a list of scopes for which this user has access.

**This endpoint returns a list of all scopes that this user has access to.**

API Keys can be used to authenticate the use of [SendGrids v3 Web API](https://sendgrid.com/docs/API_Reference/Web_API_v3/index.html), or the [Mail API Endpoint](https://sendgrid.com/docs/API_Reference/Web_API/mail.html). API Keys may be assigned certain permissions, or scopes, that limit which API endpoints they are able to access. For a more detailed explanation of how you can use API Key permissios, please visit our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/api_keys.html#-API-Key-Permissions) or [Classroom](https://sendgrid.com/docs/Classroom/Basics/API/api_key_permissions.html).

### GET /scopes


```php
$response = $sg->client->scopes()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="senders"></a>
# SENDERS

## Create a Sender Identity

**This endpoint allows you to create a new sender identity.**

*You may create up to 100 unique sender identities.*

Sender Identities are required to be verified before use. If your domain has been whitelabeled it will auto verify on creation. Otherwise an email will be sent to the `from.email`.

### POST /senders


```php
$request_body = json_decode('{
  "address": "123 Elm St.",
  "address_2": "Apt. 456",
  "city": "Denver",
  "country": "United States",
  "from": {
    "email": "from@example.com",
    "name": "Example INC"
  },
  "nickname": "My Sender ID",
  "reply_to": {
    "email": "replyto@example.com",
    "name": "Example INC"
  },
  "state": "Colorado",
  "zip": "80202"
}');
$response = $sg->client->senders()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Get all Sender Identities

**This endpoint allows you to retrieve a list of all sender identities that have been created for your account.**

Sender Identities are required to be verified before use. If your domain has been whitelabeled it will auto verify on creation. Otherwise an email will be sent to the `from.email`.

### GET /senders


```php
$response = $sg->client->senders()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a Sender Identity

**This endpoint allows you to update a sender identity.**

Updates to `from.email` require re-verification. If your domain has been whitelabeled it will auto verify on creation. Otherwise an email will be sent to the `from.email`.

Partial updates are allowed, but fields that are marked as "required" in the POST (create) endpoint must not be nil if that field is included in the PATCH request.

### PATCH /senders/{sender_id}


```php
$request_body = json_decode('{
  "address": "123 Elm St.",
  "address_2": "Apt. 456",
  "city": "Denver",
  "country": "United States",
  "from": {
    "email": "from@example.com",
    "name": "Example INC"
  },
  "nickname": "My Sender ID",
  "reply_to": {
    "email": "replyto@example.com",
    "name": "Example INC"
  },
  "state": "Colorado",
  "zip": "80202"
}');
$sender_id = "test_url_param";
$response = $sg->client->senders()->_($sender_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## View a Sender Identity

**This endpoint allows you to retrieve a specific sender identity.**

Sender Identities are required to be verified before use. If your domain has been whitelabeled it will auto verify on creation. Otherwise an email will be sent to the `from.email`.

### GET /senders/{sender_id}


```php
$sender_id = "test_url_param";
$response = $sg->client->senders()->_($sender_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a Sender Identity

**This endoint allows you to delete one of your sender identities.**

Sender Identities are required to be verified before use. If your domain has been whitelabeled it will auto verify on creation. Otherwise an email will be sent to the `from.email`.

### DELETE /senders/{sender_id}


```php
$sender_id = "test_url_param";
$response = $sg->client->senders()->_($sender_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Resend Sender Identity Verification

**This enpdoint allows you to resend a sender identity verification email.**

Sender Identities are required to be verified before use. If your domain has been whitelabeled it will auto verify on creation. Otherwise an email will be sent to the `from.email`.

### POST /senders/{sender_id}/resend_verification


```php
$sender_id = "test_url_param";
$response = $sg->client->senders()->_($sender_id)->resend_verification()->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="stats"></a>
# STATS

## Retrieve global email statistics

**This endpoint allows you to retrieve all of your global email statistics between a given date range.**

Parent accounts will see aggregated stats for their account and all subuser accounts. Subuser accounts will only see their own stats.

### GET /stats


```php
$query_params = json_decode('{"aggregated_by": "day", "limit": 1, "start_date": "2016-01-01", "end_date": "2016-04-01", "offset": 1}');
$response = $sg->client->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="subusers"></a>
# SUBUSERS

## Create Subuser

This endpoint allows you to retrieve a list of all of your subusers. You can choose to retrieve specific subusers as well as limit the results that come back from the API.

For more information about Subusers:

* [User Guide > Subusers](https://sendgrid.com/docs/User_Guide/Settings/Subusers/index.html)
* [Classroom > How do I add more subusers to my account?](https://sendgrid.com/docs/Classroom/Basics/Account/how_do_i_add_more_subusers_to_my_account.html)

### POST /subusers


```php
$request_body = json_decode('{
  "email": "John@example.com",
  "ips": [
    "1.1.1.1",
    "2.2.2.2"
  ],
  "password": "johns_password",
  "username": "John@example.com"
}');
$response = $sg->client->subusers()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## List all Subusers

This endpoint allows you to retrieve a list of all of your subusers. You can choose to retrieve specific subusers as well as limit the results that come back from the API.

For more information about Subusers:

* [User Guide > Subusers](https://sendgrid.com/docs/User_Guide/Settings/Subusers/index.html)
* [Classroom > How do I add more subusers to my account?](https://sendgrid.com/docs/Classroom/Basics/Account/how_do_i_add_more_subusers_to_my_account.html)

### GET /subusers


```php
$query_params = json_decode('{"username": "test_string", "limit": 1, "offset": 1}');
$response = $sg->client->subusers()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve Subuser Reputations

Subuser sender reputations give a good idea how well a sender is doing with regards to how recipients and recipient servers react to the mail that is being received. When a bounce, spam report, or other negative action happens on a sent email, it will effect your sender rating.

This endpoint allows you to request the reputations for your subusers.

### GET /subusers/reputations


```php
$query_params = json_decode('{"usernames": "test_string"}');
$response = $sg->client->subusers()->reputations()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve email statistics for your subusers.

**This endpoint allows you to retrieve the email statistics for the given subusers.**

You may retrieve statistics for up to 10 different subusers by including an additional _subusers_ parameter for each additional subuser.

While you can always view the statistics for all email activity on your account, subuser statistics enable you to view specific segments of your stats. Emails sent, bounces, and spam reports are always tracked for subusers. Unsubscribes, clicks, and opens are tracked if you have enabled the required settings.

For more information, see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/subuser.html).

### GET /subusers/stats


```php
$query_params = json_decode('{"end_date": "2016-04-01", "aggregated_by": "day", "limit": 1, "offset": 1, "start_date": "2016-01-01", "subusers": "test_string"}');
$response = $sg->client->subusers()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve monthly stats for all subusers

**This endpoint allows you to retrieve the monthly email statistics for all subusers over the given date range.**

While you can always view the statistics for all email activity on your account, subuser statistics enable you to view specific segments of your stats for your subusers. Emails sent, bounces, and spam reports are always tracked for subusers. Unsubscribes, clicks, and opens are tracked if you have enabled the required settings.

When using the `sort_by_metric` to sort your stats by a specific metric, you can not sort by the following metrics:
`bounce_drops`, `deferred`, `invalid_emails`, `processed`, `spam_report_drops`, `spam_reports`, or `unsubscribe_drops`.

For more information, see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/subuser.html).

### GET /subusers/stats/monthly


```php
$query_params = json_decode('{"subuser": "test_string", "limit": 1, "sort_by_metric": "test_string", "offset": 1, "date": "test_string", "sort_by_direction": "asc"}');
$response = $sg->client->subusers()->stats()->monthly()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
##  Retrieve the totals for each email statistic metric for all subusers.

**This endpoint allows you to retrieve the total sums of each email statistic metric for all subusers over the given date range.**


While you can always view the statistics for all email activity on your account, subuser statistics enable you to view specific segments of your stats. Emails sent, bounces, and spam reports are always tracked for subusers. Unsubscribes, clicks, and opens are tracked if you have enabled the required settings.

For more information, see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/subuser.html).

### GET /subusers/stats/sums


```php
$query_params = json_decode('{"end_date": "2016-04-01", "aggregated_by": "day", "limit": 1, "sort_by_metric": "test_string", "offset": 1, "start_date": "2016-01-01", "sort_by_direction": "asc"}');
$response = $sg->client->subusers()->stats()->sums()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Enable/disable a subuser

This endpoint allows you to enable or disable a subuser.

For more information about Subusers:

* [User Guide > Subusers](https://sendgrid.com/docs/User_Guide/Settings/Subusers/index.html)
* [Classroom > How do I add more subusers to my account?](https://sendgrid.com/docs/Classroom/Basics/Account/how_do_i_add_more_subusers_to_my_account.html)

### PATCH /subusers/{subuser_name}


```php
$request_body = json_decode('{
  "disabled": false
}');
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a subuser

This endpoint allows you to delete a subuser. This is a permanent action, once deleted a subuser cannot be retrieved.

For more information about Subusers:

* [User Guide > Subusers](https://sendgrid.com/docs/User_Guide/Settings/Subusers/index.html)
* [Classroom > How do I add more subusers to my account?](https://sendgrid.com/docs/Classroom/Basics/Account/how_do_i_add_more_subusers_to_my_account.html)

### DELETE /subusers/{subuser_name}


```php
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update IPs assigned to a subuser

Each subuser should be assigned to an IP address, from which all of this subuser's mail will be sent. Often, this is the same IP as the parent account, but each subuser can have their own, or multiple, IP addresses as well.

More information:

* [How to request more IPs](https://sendgrid.com/docs/Classroom/Basics/Account/adding_an_additional_dedicated_ip_to_your_account.html)
* [IPs can be whitelabeled](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/ips.html)

### PUT /subusers/{subuser_name}/ips


```php
$request_body = json_decode('[
  "127.0.0.1"
]');
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->ips()->put($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Monitor Settings for a subuser

Subuser monitor settings allow you to receive a sample of an outgoing message by a specific customer at a specific frequency of emails.

### PUT /subusers/{subuser_name}/monitor


```php
$request_body = json_decode('{
  "email": "example@example.com",
  "frequency": 500
}');
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->monitor()->put($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create monitor settings

Subuser monitor settings allow you to receive a sample of an outgoing message by a specific customer at a specific frequency of emails.

### POST /subusers/{subuser_name}/monitor


```php
$request_body = json_decode('{
  "email": "example@example.com",
  "frequency": 50000
}');
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->monitor()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve monitor settings for a subuser

Subuser monitor settings allow you to receive a sample of an outgoing message by a specific customer at a specific frequency of emails.

### GET /subusers/{subuser_name}/monitor


```php
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->monitor()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete monitor settings

Subuser monitor settings allow you to receive a sample of an outgoing message by a specific customer at a specific frequency of emails.

### DELETE /subusers/{subuser_name}/monitor


```php
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->monitor()->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve the monthly email statistics for a single subuser

**This endpoint allows you to retrive the monthly email statistics for a specific subuser.**

While you can always view the statistics for all email activity on your account, subuser statistics enable you to view specific segments of your stats for your subusers. Emails sent, bounces, and spam reports are always tracked for subusers. Unsubscribes, clicks, and opens are tracked if you have enabled the required settings.

When using the `sort_by_metric` to sort your stats by a specific metric, you can not sort by the following metrics:
`bounce_drops`, `deferred`, `invalid_emails`, `processed`, `spam_report_drops`, `spam_reports`, or `unsubscribe_drops`.

For more information, see our [User Guide](https://sendgrid.com/docs/User_Guide/Statistics/subuser.html).

### GET /subusers/{subuser_name}/stats/monthly


```php
$query_params = json_decode('{"date": "test_string", "sort_by_direction": "asc", "limit": 1, "sort_by_metric": "test_string", "offset": 1}');
$subuser_name = "test_url_param";
$response = $sg->client->subusers()->_($subuser_name)->stats()->monthly()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="suppression"></a>
# SUPPRESSION

## Retrieve all blocks

**This endpoint allows you to retrieve a list of all email addresses that are currently on your blocks list.**

[Blocks](https://sendgrid.com/docs/Glossary/blocks.html) happen when your message was rejected for a reason related to the message, not the recipient address. This can happen when your mail server IP address has been added to a blacklist or blocked by an ISP, or if the message content is flagged by a filter on the receiving server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/blocks.html).

### GET /suppression/blocks


```php
$query_params = json_decode('{"start_time": 1, "limit": 1, "end_time": 1, "offset": 1}');
$response = $sg->client->suppression()->blocks()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete blocks

**This endpoint allows you to delete all email addresses on your blocks list.**

There are two options for deleting blocked emails:

1. You can delete all blocked emails by setting `delete_all` to true in the request body.
2. You can delete some blocked emails by specifying the email addresses in an array in the request body.

[Blocks](https://sendgrid.com/docs/Glossary/blocks.html) happen when your message was rejected for a reason related to the message, not the recipient address. This can happen when your mail server IP address has been added to a blacklist or blocked by an ISP, or if the message content is flagged by a filter on the receiving server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/blocks.html).

### DELETE /suppression/blocks


```php
$request_body = json_decode('{
  "delete_all": false,
  "emails": [
    "example1@example.com",
    "example2@example.com"
  ]
}');
$response = $sg->client->suppression()->blocks()->delete($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a specific block

**This endpoint allows you to retrieve a specific email address from your blocks list.**

[Blocks](https://sendgrid.com/docs/Glossary/blocks.html) happen when your message was rejected for a reason related to the message, not the recipient address. This can happen when your mail server IP address has been added to a blacklist or blocked by an ISP, or if the message content is flagged by a filter on the receiving server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/blocks.html).

### GET /suppression/blocks/{email}


```php
$email = "test_url_param";
$response = $sg->client->suppression()->blocks()->_($email)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a specific block

**This endpoint allows you to delete a specific email address from your blocks list.**

[Blocks](https://sendgrid.com/docs/Glossary/blocks.html) happen when your message was rejected for a reason related to the message, not the recipient address. This can happen when your mail server IP address has been added to a blacklist or blocked by an ISP, or if the message content is flagged by a filter on the receiving server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/blocks.html).

### DELETE /suppression/blocks/{email}


```php
$email = "test_url_param";
$response = $sg->client->suppression()->blocks()->_($email)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all bounces

**This endpoint allows you to retrieve all of your bounces.**

Bounces are messages that are returned to the server that sent it.

For more information see:

* [User Guide > Bounces](https://sendgrid.com/docs/User_Guide/Suppressions/bounces.html) for more information
* [Glossary > Bounces](https://sendgrid.com/docs/Glossary/Bounces.html)

### GET /suppression/bounces


```php
$query_params = json_decode('{"start_time": 1, "end_time": 1}');
$response = $sg->client->suppression()->bounces()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete bounces

**This endpoint allows you to delete all of your bounces. You can also use this endpoint to remove a specific email address from your bounce list.**

Bounces are messages that are returned to the server that sent it.

For more information see:

* [User Guide > Bounces](https://sendgrid.com/docs/User_Guide/Suppressions/bounces.html) for more information
* [Glossary > Bounces](https://sendgrid.com/docs/Glossary/Bounces.html)
* [Classroom > List Scrubbing Guide](https://sendgrid.com/docs/Classroom/Deliver/list_scrubbing.html)

Note: the `delete_all` and `emails` parameters should be used independently of each other as they have different purposes.

### DELETE /suppression/bounces


```php
$request_body = json_decode('{
  "delete_all": true,
  "emails": [
    "example@example.com",
    "example2@example.com"
  ]
}');
$response = $sg->client->suppression()->bounces()->delete($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a Bounce

**This endpoint allows you to retrieve a specific bounce for a given email address.**

Bounces are messages that are returned to the server that sent it.

For more information see:

* [User Guide > Bounces](https://sendgrid.com/docs/User_Guide/Suppressions/bounces.html) for more information
* [Glossary > Bounces](https://sendgrid.com/docs/Glossary/Bounces.html)
* [Classroom > List Scrubbing Guide](https://sendgrid.com/docs/Classroom/Deliver/list_scrubbing.html)

### GET /suppression/bounces/{email}


```php
$email = "test_url_param";
$response = $sg->client->suppression()->bounces()->_($email)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a bounce

**This endpoint allows you to remove an email address from your bounce list.**

Bounces are messages that are returned to the server that sent it. This endpoint allows you to delete a single email addresses from your bounce list.

For more information see:

* [User Guide > Bounces](https://sendgrid.com/docs/User_Guide/Suppressions/bounces.html) for more information
* [Glossary > Bounces](https://sendgrid.com/docs/Glossary/Bounces.html)
* [Classroom > List Scrubbing Guide](https://sendgrid.com/docs/Classroom/Deliver/list_scrubbing.html)

### DELETE /suppression/bounces/{email}


```php
$query_params = json_decode('{"email_address": "example@example.com"}');
$email = "test_url_param";
$response = $sg->client->suppression()->bounces()->_($email)->delete(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all invalid emails

**This endpoint allows you to retrieve a list of all invalid email addresses.**

An invalid email occurs when you attempt to send email to an address that is formatted in a manner that does not meet internet email format standards or the email does not exist at the recipients mail server.

Examples include addresses without the @ sign or addresses that include certain special characters and/or spaces. This response can come from our own server or the recipient mail server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/invalid_emails.html).

### GET /suppression/invalid_emails


```php
$query_params = json_decode('{"start_time": 1, "limit": 1, "end_time": 1, "offset": 1}');
$response = $sg->client->suppression()->invalid_emails()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete invalid emails

**This endpoint allows you to remove email addresses from your invalid email address list.**

There are two options for deleting invalid email addresses:

1) You can delete all invalid email addresses by setting `delete_all` to true in the request body.
2) You can delete some invalid email addresses by specifying certain addresses in an array in the request body.

An invalid email occurs when you attempt to send email to an address that is formatted in a manner that does not meet internet email format standards or the email does not exist at the recipients mail server.

Examples include addresses without the @ sign or addresses that include certain special characters and/or spaces. This response can come from our own server or the recipient mail server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/invalid_emails.html).

### DELETE /suppression/invalid_emails


```php
$request_body = json_decode('{
  "delete_all": false,
  "emails": [
    "example1@example.com",
    "example2@example.com"
  ]
}');
$response = $sg->client->suppression()->invalid_emails()->delete($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a specific invalid email

**This endpoint allows you to retrieve a specific invalid email addresses.**

An invalid email occurs when you attempt to send email to an address that is formatted in a manner that does not meet internet email format standards or the email does not exist at the recipients mail server.

Examples include addresses without the @ sign or addresses that include certain special characters and/or spaces. This response can come from our own server or the recipient mail server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/invalid_emails.html).

### GET /suppression/invalid_emails/{email}


```php
$email = "test_url_param";
$response = $sg->client->suppression()->invalid_emails()->_($email)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a specific invalid email

**This endpoint allows you to remove a specific email address from the invalid email address list.**

An invalid email occurs when you attempt to send email to an address that is formatted in a manner that does not meet internet email format standards or the email does not exist at the recipients mail server.

Examples include addresses without the @ sign or addresses that include certain special characters and/or spaces. This response can come from our own server or the recipient mail server.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/invalid_emails.html).

### DELETE /suppression/invalid_emails/{email}


```php
$email = "test_url_param";
$response = $sg->client->suppression()->invalid_emails()->_($email)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a specific spam report

**This endpoint allows you to retrieve a specific spam report.**

[Spam reports](https://sendgrid.com/docs/Glossary/spam_reports.html) happen when a recipient indicates that they think your email is [spam](https://sendgrid.com/docs/Glossary/spam.html) and then their email provider reports this to SendGrid.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/spam_reports.html).

### GET /suppression/spam_report/{email}


```php
$email = "test_url_param";
$response = $sg->client->suppression()->spam_report()->_($email)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a specific spam report

**This endpoint allows you to delete a specific spam report.**

[Spam reports](https://sendgrid.com/docs/Glossary/spam_reports.html) happen when a recipient indicates that they think your email is [spam](https://sendgrid.com/docs/Glossary/spam.html) and then their email provider reports this to SendGrid.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/spam_reports.html).

### DELETE /suppression/spam_report/{email}


```php
$email = "test_url_param";
$response = $sg->client->suppression()->spam_report()->_($email)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all spam reports

**This endpoint allows you to retrieve all spam reports.**

[Spam reports](https://sendgrid.com/docs/Glossary/spam_reports.html) happen when a recipient indicates that they think your email is [spam](https://sendgrid.com/docs/Glossary/spam.html) and then their email provider reports this to SendGrid.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/spam_reports.html).

### GET /suppression/spam_reports


```php
$query_params = json_decode('{"start_time": 1, "limit": 1, "end_time": 1, "offset": 1}');
$response = $sg->client->suppression()->spam_reports()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete spam reports

**This endpoint allows you to delete your spam reports.**

There are two options for deleting spam reports:

1) You can delete all spam reports by setting "delete_all" to true in the request body.
2) You can delete some spam reports by specifying the email addresses in an array in the request body.

[Spam reports](https://sendgrid.com/docs/Glossary/spam_reports.html) happen when a recipient indicates that they think your email is [spam](https://sendgrid.com/docs/Glossary/spam.html) and then their email provider reports this to SendGrid.

For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/spam_reports.html).

### DELETE /suppression/spam_reports


```php
$request_body = json_decode('{
  "delete_all": false,
  "emails": [
    "example1@example.com",
    "example2@example.com"
  ]
}');
$response = $sg->client->suppression()->spam_reports()->delete($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all global suppressions

**This endpoint allows you to retrieve a list of all email address that are globally suppressed.**

A global suppression (or global unsubscribe) is an email address of a recipient who does not want to receive any of your messages. A globally suppressed recipient will be removed from any email you send. For more information, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Suppressions/global_unsubscribes.html).

### GET /suppression/unsubscribes


```php
$query_params = json_decode('{"start_time": 1, "limit": 1, "end_time": 1, "offset": 1}');
$response = $sg->client->suppression()->unsubscribes()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="templates"></a>
# TEMPLATES

## Create a transactional template.

**This endpoint allows you to create a transactional template.**

Each user can create up to 300 different transactional templates. Transactional templates are specific to accounts and subusers. Templates created on a parent account will not be accessible from the subuser accounts.

Transactional templates are templates created specifically for transactional email and are not to be confused with [Marketing Campaigns templates](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/templates.html). For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

### POST /templates


```php
$request_body = json_decode('{
  "name": "example_name"
}');
$response = $sg->client->templates()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all transactional templates.

**This endpoint allows you to retrieve all transactional templates.**

Each user can create up to 300 different transactional templates. Transactional templates are specific to accounts and subusers. Templates created on a parent account will not be accessible from the subuser accounts.

Transactional templates are templates created specifically for transactional email and are not to be confused with [Marketing Campaigns templates](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/templates.html). For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

### GET /templates


```php
$response = $sg->client->templates()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Edit a transactional template.

**This endpoint allows you to edit a transactional template.**

Each user can create up to 300 different transactional templates. Transactional templates are specific to accounts and subusers. Templates created on a parent account will not be accessible from the subuser accounts.

Transactional templates are templates created specifically for transactional email and are not to be confused with [Marketing Campaigns templates](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/templates.html). For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).


### PATCH /templates/{template_id}


```php
$request_body = json_decode('{
  "name": "new_example_name"
}');
$template_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a single transactional template.

**This endpoint allows you to retrieve a single transactional template.**

Each user can create up to 300 different transactional templates. Transactional templates are specific to accounts and subusers. Templates created on a parent account will not be accessible from the subuser accounts.

Transactional templates are templates created specifically for transactional email and are not to be confused with [Marketing Campaigns templates](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/templates.html). For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).


### GET /templates/{template_id}


```php
$template_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a template.

**This endpoint allows you to delete a transactional template.**

Each user can create up to 300 different transactional templates. Transactional templates are specific to accounts and subusers. Templates created on a parent account will not be accessible from the subuser accounts.

Transactional templates are templates created specifically for transactional email and are not to be confused with [Marketing Campaigns templates](https://sendgrid.com/docs/User_Guide/Marketing_Campaigns/templates.html). For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).


### DELETE /templates/{template_id}


```php
$template_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create a new transactional template version.

**This endpoint allows you to create a new version of a template.**

Each transactional template can have multiple versions, each version with its own subject and content. Each user can have up to 300 versions across across all templates.

For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).


### POST /templates/{template_id}/versions


```php
$request_body = json_decode('{
  "active": 1,
  "html_content": "<%body%>",
  "name": "example_version_name",
  "plain_content": "<%body%>",
  "subject": "<%subject%>",
  "template_id": "ddb96bbc-9b92-425e-8979-99464621b543"
}');
$template_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->versions()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Edit a transactional template version.

**This endpoint allows you to edit a version of one of your transactional templates.**

Each transactional template can have multiple versions, each version with its own subject and content. Each user can have up to 300 versions across across all templates.

For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

## URI Parameters
| URI Parameter | Type | Description |
|---|---|---|
| template_id | string | The ID of the original template |
| version_id | string | The ID of the template version |

### PATCH /templates/{template_id}/versions/{version_id}


```php
$request_body = json_decode('{
  "active": 1,
  "html_content": "<%body%>",
  "name": "updated_example_name",
  "plain_content": "<%body%>",
  "subject": "<%subject%>"
}');
$template_id = "test_url_param";
$version_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->versions()->_($version_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a specific transactional template version.

**This endpoint allows you to retrieve a specific version of a template.**

Each transactional template can have multiple versions, each version with its own subject and content. Each user can have up to 300 versions across across all templates.

For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

## URI Parameters
| URI Parameter | Type | Description |
|---|---|---|
| template_id | string | The ID of the original template |
| version_id | string |  The ID of the template version |

### GET /templates/{template_id}/versions/{version_id}


```php
$template_id = "test_url_param";
$version_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->versions()->_($version_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a transactional template version.

**This endpoint allows you to delete one of your transactional template versions.**

Each transactional template can have multiple versions, each version with its own subject and content. Each user can have up to 300 versions across across all templates.

For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

## URI Parameters
| URI Parameter | Type | Description |
|---|---|---|
| template_id | string | The ID of the original template |
| version_id | string | The ID of the template version |

### DELETE /templates/{template_id}/versions/{version_id}


```php
$template_id = "test_url_param";
$version_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->versions()->_($version_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Activate a transactional template version.

**This endpoint allows you to activate a version of one of your templates.**

Each transactional template can have multiple versions, each version with its own subject and content. Each user can have up to 300 versions across across all templates.


For more information about transactional templates, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html).

## URI Parameters
| URI Parameter | Type | Description |
|---|---|---|
| template_id | string | The ID of the original template |
| version_id | string |  The ID of the template version |

### POST /templates/{template_id}/versions/{version_id}/activate


```php
$template_id = "test_url_param";
$version_id = "test_url_param";
$response = $sg->client->templates()->_($template_id)->versions()->_($version_id)->activate()->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="tracking_settings"></a>
# TRACKING SETTINGS

## Retrieve Tracking Settings

**This endpoint allows you to retrieve a list of all tracking settings that you can enable on your account.**

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### GET /tracking_settings


```php
$query_params = json_decode('{"limit": 1, "offset": 1}');
$response = $sg->client->tracking_settings()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Click Tracking Settings

**This endpoint allows you to change your current click tracking setting. You can enable, or disable, click tracking using this endpoint.**

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### PATCH /tracking_settings/click


```php
$request_body = json_decode('{
  "enabled": true
}');
$response = $sg->client->tracking_settings()->click()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve Click Track Settings

**This endpoint allows you to retrieve your current click tracking setting.**

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### GET /tracking_settings/click


```php
$response = $sg->client->tracking_settings()->click()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Google Analytics Settings

**This endpoint allows you to update your current setting for Google Analytics.**

For more information about using Google Analytics, please refer to [Googles URL Builder](https://support.google.com/analytics/answer/1033867?hl=en) and their article on ["Best Practices for Campaign Building"](https://support.google.com/analytics/answer/1037445).

We default the settings to Googles recommendations. For more information, see [Google Analytics Demystified](https://sendgrid.com/docs/Classroom/Track/Collecting_Data/google_analytics_demystified_ga_statistics_vs_sg_statistics.html).

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### PATCH /tracking_settings/google_analytics


```php
$request_body = json_decode('{
  "enabled": true,
  "utm_campaign": "website",
  "utm_content": "",
  "utm_medium": "email",
  "utm_source": "sendgrid.com",
  "utm_term": ""
}');
$response = $sg->client->tracking_settings()->google_analytics()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve Google Analytics Settings

**This endpoint allows you to retrieve your current setting for Google Analytics.**

For more information about using Google Analytics, please refer to [Googles URL Builder](https://support.google.com/analytics/answer/1033867?hl=en) and their article on ["Best Practices for Campaign Building"](https://support.google.com/analytics/answer/1037445).

We default the settings to Googles recommendations. For more information, see [Google Analytics Demystified](https://sendgrid.com/docs/Classroom/Track/Collecting_Data/google_analytics_demystified_ga_statistics_vs_sg_statistics.html).

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### GET /tracking_settings/google_analytics


```php
$response = $sg->client->tracking_settings()->google_analytics()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Open Tracking Settings

**This endpoint allows you to update your current settings for open tracking.**

Open Tracking adds an invisible image at the end of the email which can track email opens. If the email recipient has images enabled on their email client, a request to SendGrids server for the invisible image is executed and an open event is logged. These events are logged in the Statistics portal, Email Activity interface, and are reported by the Event Webhook.

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### PATCH /tracking_settings/open


```php
$request_body = json_decode('{
  "enabled": true
}');
$response = $sg->client->tracking_settings()->open()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Get Open Tracking Settings

**This endpoint allows you to retrieve your current settings for open tracking.**

Open Tracking adds an invisible image at the end of the email which can track email opens. If the email recipient has images enabled on their email client, a request to SendGrids server for the invisible image is executed and an open event is logged. These events are logged in the Statistics portal, Email Activity interface, and are reported by the Event Webhook.

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### GET /tracking_settings/open


```php
$response = $sg->client->tracking_settings()->open()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Subscription Tracking Settings

**This endpoint allows you to update your current settings for subscription tracking.**

Subscription tracking adds links to the bottom of your emails that allows your recipients to subscribe to, or unsubscribe from, your emails.

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### PATCH /tracking_settings/subscription


```php
$request_body = json_decode('{
  "enabled": true,
  "html_content": "html content",
  "landing": "landing page html",
  "plain_content": "text content",
  "replace": "replacement tag",
  "url": "url"
}');
$response = $sg->client->tracking_settings()->subscription()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve Subscription Tracking Settings

**This endpoint allows you to retrieve your current settings for subscription tracking.**

Subscription tracking adds links to the bottom of your emails that allows your recipients to subscribe to, or unsubscribe from, your emails.

You can track a variety of the actions your recipients may take when interacting with your emails including opening your emails, clicking on links in your emails, and subscribing to (or unsubscribing from) your emails.

For more information about tracking, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/tracking.html).

### GET /tracking_settings/subscription


```php
$response = $sg->client->tracking_settings()->subscription()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="user"></a>
# USER

## Get a user's account information.

**This endpoint allows you to retrieve your user account details.**

Your user's account information includes the user's account type and reputation.

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

### GET /user/account


```php
$response = $sg->client->user()->account()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve your credit balance

**This endpoint allows you to retrieve the current credit balance for your account.**

Your monthly credit allotment limits the number of emails you may send before incurring overage charges. For more information about credits and billing, please visit our [Clssroom](https://sendgrid.com/docs/Classroom/Basics/Billing/billing_info_and_faqs.html).

### GET /user/credits


```php
$response = $sg->client->user()->credits()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update your account email address

**This endpoint allows you to update the email address currently on file for your account.**

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

### PUT /user/email


```php
$request_body = json_decode('{
  "email": "example@example.com"
}');
$response = $sg->client->user()->email()->put($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve your account email address

**This endpoint allows you to retrieve the email address currently on file for your account.**

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

### GET /user/email


```php
$response = $sg->client->user()->email()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update your password

**This endpoint allows you to update your password.**

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

### PUT /user/password


```php
$request_body = json_decode('{
  "new_password": "new_password",
  "old_password": "old_password"
}');
$response = $sg->client->user()->password()->put($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a user's profile

**This endpoint allows you to update your current profile details.**

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

It should be noted that any one or more of the parameters can be updated via the PATCH /user/profile endpoint. The only requirement is that you include at least one when you PATCH.

### PATCH /user/profile


```php
$request_body = json_decode('{
  "city": "Orange",
  "first_name": "Example",
  "last_name": "User"
}');
$response = $sg->client->user()->profile()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Get a user's profile

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

### GET /user/profile


```php
$response = $sg->client->user()->profile()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Cancel or pause a scheduled send

**This endpoint allows you to cancel or pause an email that has been scheduled to be sent.**

If the maximum number of cancellations/pauses are added, HTTP 400 will
be returned.

The Cancel Scheduled Sends feature allows the customer to cancel a scheduled send based on a Batch ID included in the SMTPAPI header.Scheduled sends cancelled less than 10 minutes before the scheduled time are not guaranteed to be cancelled.

### POST /user/scheduled_sends


```php
$request_body = json_decode('{
  "batch_id": "YOUR_BATCH_ID",
  "status": "pause"
}');
$response = $sg->client->user()->scheduled_sends()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all scheduled sends

**This endpoint allows you to retrieve all cancel/paused scheduled send information.**

The Cancel Scheduled Sends feature allows the customer to cancel a scheduled send based on a Batch ID included in the SMTPAPI header.Scheduled sends cancelled less than 10 minutes before the scheduled time are not guaranteed to be cancelled.

### GET /user/scheduled_sends


```php
$response = $sg->client->user()->scheduled_sends()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update user scheduled send information

**This endpoint allows you to update the status of a scheduled send for the given `batch_id`.**

The Cancel Scheduled Sends feature allows the customer to cancel a scheduled send based on a Batch ID included in the SMTPAPI header.Scheduled sends cancelled less than 10 minutes before the scheduled time are not guaranteed to be cancelled.

### PATCH /user/scheduled_sends/{batch_id}


```php
$request_body = json_decode('{
  "status": "pause"
}');
$batch_id = "test_url_param";
$response = $sg->client->user()->scheduled_sends()->_($batch_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve scheduled send

**This endpoint allows you to retrieve the cancel/paused scheduled send information for a specific `batch_id`.**

The Cancel Scheduled Sends feature allows the customer to cancel a scheduled send based on a Batch ID included in the SMTPAPI header.Scheduled sends cancelled less than 10 minutes before the scheduled time are not guaranteed to be cancelled.

### GET /user/scheduled_sends/{batch_id}


```php
$batch_id = "test_url_param";
$response = $sg->client->user()->scheduled_sends()->_($batch_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a cancellation or pause of a scheduled send

**This endpoint allows you to delete the cancellation/pause of a scheduled send.**

The Cancel Scheduled Sends feature allows the customer to cancel a scheduled send based on a Batch ID included in the SMTPAPI header.Scheduled sends cancelled less than 10 minutes before the scheduled time are not guaranteed to be cancelled.

### DELETE /user/scheduled_sends/{batch_id}


```php
$batch_id = "test_url_param";
$response = $sg->client->user()->scheduled_sends()->_($batch_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Enforced TLS settings

**This endpoint allows you to update your current Enforced TLS settings.**

The Enforced TLS settings specify whether or not the recipient is required to support TLS or have a valid certificate. See the [SMTP Ports User Guide](https://sendgrid.com/docs/Classroom/Basics/Email_Infrastructure/smtp_ports.html) for more information on opportunistic TLS.

**Note:** If either setting is enabled and the recipient does not support TLS or have a valid certificate, we drop the message and send a block event with TLS required but not supported as the description.

### PATCH /user/settings/enforced_tls


```php
$request_body = json_decode('{
  "require_tls": true,
  "require_valid_cert": false
}');
$response = $sg->client->user()->settings()->enforced_tls()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve current Enforced TLS settings.

**This endpoint allows you to retrieve your current Enforced TLS settings.**

The Enforced TLS settings specify whether or not the recipient is required to support TLS or have a valid certificate. See the [SMTP Ports User Guide](https://sendgrid.com/docs/Classroom/Basics/Email_Infrastructure/smtp_ports.html) for more information on opportunistic TLS.

**Note:** If either setting is enabled and the recipient does not support TLS or have a valid certificate, we drop the message and send a block event with TLS required but not supported as the description.

### GET /user/settings/enforced_tls


```php
$response = $sg->client->user()->settings()->enforced_tls()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update your username

**This endpoint allows you to update the username for your account.**

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

### PUT /user/username


```php
$request_body = json_decode('{
  "username": "test_username"
}');
$response = $sg->client->user()->username()->put($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve your username

**This endpoint allows you to retrieve your current account username.**

Keeping your user profile up to date is important. This will help SendGrid to verify who you are as well as contact you should we need to.

For more information about your user profile:

* [SendGrid Account Settings](https://sendgrid.com/docs/User_Guide/Settings/account.html)

### GET /user/username


```php
$response = $sg->client->user()->username()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update Event Notification Settings

**This endpoint allows you to update your current event webhook settings.**

If an event type is marked as `true`, then the event webhook will include information about that event.

SendGrids Event Webhook will notify a URL of your choice via HTTP POST with information about events that occur as SendGrid processes your email.

Common uses of this data are to remove unsubscribes, react to spam reports, determine unengaged recipients, identify bounced email addresses, or create advanced analytics of your email program.

### PATCH /user/webhooks/event/settings


```php
$request_body = json_decode('{
  "bounce": true,
  "click": true,
  "deferred": true,
  "delivered": true,
  "dropped": true,
  "enabled": true,
  "group_resubscribe": true,
  "group_unsubscribe": true,
  "open": true,
  "processed": true,
  "spam_report": true,
  "unsubscribe": true,
  "url": "url"
}');
$response = $sg->client->user()->webhooks()->event()->settings()->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve Event Webhook settings

**This endpoint allows you to retrieve your current event webhook settings.**

If an event type is marked as `true`, then the event webhook will include information about that event.

SendGrids Event Webhook will notify a URL of your choice via HTTP POST with information about events that occur as SendGrid processes your email.

Common uses of this data are to remove unsubscribes, react to spam reports, determine unengaged recipients, identify bounced email addresses, or create advanced analytics of your email program.

### GET /user/webhooks/event/settings


```php
$response = $sg->client->user()->webhooks()->event()->settings()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Test Event Notification Settings

**This endpoint allows you to test your event webhook by sending a fake event notification post to the provided URL.**

SendGrids Event Webhook will notify a URL of your choice via HTTP POST with information about events that occur as SendGrid processes your email.

Common uses of this data are to remove unsubscribes, react to spam reports, determine unengaged recipients, identify bounced email addresses, or create advanced analytics of your email program.

### POST /user/webhooks/event/test


```php
$request_body = json_decode('{
  "url": "url"
}');
$response = $sg->client->user()->webhooks()->event()->test()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create a parse setting

**This endpoint allows you to create a new inbound parse setting.**

The inbound parse webhook allows you to have incoming emails parsed, extracting some or all of the content, and then have that content POSTed by SendGrid to a URL of your choosing. For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Webhooks/parse.html).

### POST /user/webhooks/parse/settings


```php
$request_body = json_decode('{
  "hostname": "myhostname.com",
  "send_raw": false,
  "spam_check": true,
  "url": "http://email.myhosthame.com"
}');
$response = $sg->client->user()->webhooks()->parse()->settings()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all parse settings

**This endpoint allows you to retrieve all of your current inbound parse settings.**

The inbound parse webhook allows you to have incoming emails parsed, extracting some or all of the contnet, and then have that content POSTed by SendGrid to a URL of your choosing. For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Webhooks/parse.html).

### GET /user/webhooks/parse/settings


```php
$response = $sg->client->user()->webhooks()->parse()->settings()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a parse setting

**This endpoint allows you to update a specific inbound parse setting.**

The inbound parse webhook allows you to have incoming emails parsed, extracting some or all of the contnet, and then have that content POSTed by SendGrid to a URL of your choosing. For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Webhooks/parse.html).

### PATCH /user/webhooks/parse/settings/{hostname}


```php
$request_body = json_decode('{
  "send_raw": true,
  "spam_check": false,
  "url": "http://newdomain.com/parse"
}');
$hostname = "test_url_param";
$response = $sg->client->user()->webhooks()->parse()->settings()->_($hostname)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a specific parse setting

**This endpoint allows you to retrieve a specific inbound parse setting.**

The inbound parse webhook allows you to have incoming emails parsed, extracting some or all of the contnet, and then have that content POSTed by SendGrid to a URL of your choosing. For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Webhooks/parse.html).

### GET /user/webhooks/parse/settings/{hostname}


```php
$hostname = "test_url_param";
$response = $sg->client->user()->webhooks()->parse()->settings()->_($hostname)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a parse setting

**This endpoint allows you to delete a specific inbound parse setting.**

The inbound parse webhook allows you to have incoming emails parsed, extracting some or all of the contnet, and then have that content POSTed by SendGrid to a URL of your choosing. For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Webhooks/parse.html).

### DELETE /user/webhooks/parse/settings/{hostname}


```php
$hostname = "test_url_param";
$response = $sg->client->user()->webhooks()->parse()->settings()->_($hostname)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieves Inbound Parse Webhook statistics.

**This endpoint allows you to retrieve the statistics for your Parse Webhook useage.**

SendGrid's Inbound Parse Webhook allows you to parse the contents and attachments of incomming emails. The Parse API can then POST the parsed emails to a URL that you specify. The Inbound Parse Webhook cannot parse messages greater than 20MB in size, including all attachments.

There are a number of pre-made integrations for the SendGrid Parse Webhook which make processing events easy. You can find these integrations in the [Library Index](https://sendgrid.com/docs/Integrate/libraries.html#-Webhook-Libraries).

### GET /user/webhooks/parse/stats


```php
$query_params = json_decode('{"aggregated_by": "day", "limit": "test_string", "start_date": "2016-01-01", "end_date": "2016-04-01", "offset": "test_string"}');
$response = $sg->client->user()->webhooks()->parse()->stats()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
<a name="whitelabel"></a>
# WHITELABEL

## Create a domain whitelabel.

**This endpoint allows you to create a whitelabel for one of your domains.**

If you are creating a domain whitelabel that you would like a subuser to use, you have two options:
1. Use the "username" parameter. This allows you to create a whitelabel on behalf of your subuser. This means the subuser is able to see and modify the created whitelabel.
2. Use the Association workflow (see Associate Domain section). This allows you to assign a whitelabel created by the parent to a subuser. This means the subuser will default to the assigned whitelabel, but will not be able to see or modify that whitelabel. However, if the subuser creates their own whitelabel it will overwrite the assigned whitelabel.

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

### POST /whitelabel/domains


```php
$request_body = json_decode('{
  "automatic_security": false,
  "custom_spf": true,
  "default": true,
  "domain": "example.com",
  "ips": [
    "192.168.1.1",
    "192.168.1.2"
  ],
  "subdomain": "news",
  "username": "john@example.com"
}');
$response = $sg->client->whitelabel()->domains()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## List all domain whitelabels.

**This endpoint allows you to retrieve a list of all domain whitelabels you have created.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)


### GET /whitelabel/domains


```php
$query_params = json_decode('{"username": "test_string", "domain": "test_string", "exclude_subusers": "true", "limit": 1, "offset": 1}');
$response = $sg->client->whitelabel()->domains()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Get the default domain whitelabel.

**This endpoint allows you to retrieve the default whitelabel for a domain.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

## URI Parameters
| URI Parameter   | Type   | Description  |
|---|---|---|
| domain | string  |The domain to find a default domain whitelabel for. |

### GET /whitelabel/domains/default


```php
$response = $sg->client->whitelabel()->domains()->default()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## List the domain whitelabel associated with the given user.

**This endpoint allows you to retrieve all of the whitelabels that have been assigned to a specific subuser.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

Domain whitelabels can be associated with (i.e. assigned to) subusers from a parent account. This functionality allows subusers to send mail using their parent's whitelabels. To associate a whitelabel with a subuser, the parent account must first create the whitelabel and validate it. The the parent may then associate the whitelabel via the subuser management tools.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

## URI Parameters
| URI Parameter   | Type  | Description  |
|---|---|---|
| username | string  | Username of the subuser to find associated whitelabels for. |

### GET /whitelabel/domains/subuser


```php
$response = $sg->client->whitelabel()->domains()->subuser()->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Disassociate a domain whitelabel from a given user.

**This endpoint allows you to disassociate a specific whitelabel from a subuser.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

Domain whitelabels can be associated with (i.e. assigned to) subusers from a parent account. This functionality allows subusers to send mail using their parent's whitelabels. To associate a whitelabel with a subuser, the parent account must first create the whitelabel and validate it. The the parent may then associate the whitelabel via the subuser management tools.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

## URI Parameters
| URI Parameter   | Type  | Required?  | Description  |
|---|---|---|---|
| username | string  | required  | Username for the subuser to find associated whitelabels for. |

### DELETE /whitelabel/domains/subuser


```php
$response = $sg->client->whitelabel()->domains()->subuser()->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a domain whitelabel.

**This endpoint allows you to update the settings for a domain whitelabel.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

### PATCH /whitelabel/domains/{domain_id}


```php
$request_body = json_decode('{
  "custom_spf": true,
  "default": false
}');
$domain_id = "test_url_param";
$response = $sg->client->whitelabel()->domains()->_($domain_id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a domain whitelabel.

**This endpoint allows you to retrieve a specific domain whitelabel.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)


### GET /whitelabel/domains/{domain_id}


```php
$domain_id = "test_url_param";
$response = $sg->client->whitelabel()->domains()->_($domain_id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a domain whitelabel.

**This endpoint allows you to delete a domain whitelabel.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

### DELETE /whitelabel/domains/{domain_id}


```php
$domain_id = "test_url_param";
$response = $sg->client->whitelabel()->domains()->_($domain_id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Associate a domain whitelabel with a given user.

**This endpoint allows you to associate a specific domain whitelabel with a subuser.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

Domain whitelabels can be associated with (i.e. assigned to) subusers from a parent account. This functionality allows subusers to send mail using their parent's whitelabels. To associate a whitelabel with a subuser, the parent account must first create the whitelabel and validate it. The the parent may then associate the whitelabel via the subuser management tools.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

## URI Parameters
| URI Parameter   | Type   | Description  |
|---|---|---|
| domain_id | integer   | ID of the domain whitelabel to associate with the subuser. |

### POST /whitelabel/domains/{domain_id}/subuser


```php
$request_body = json_decode('{
  "username": "jane@example.com"
}');
$domain_id = "test_url_param";
$response = $sg->client->whitelabel()->domains()->_($domain_id)->subuser()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Add an IP to a domain whitelabel.

**This endpoint allows you to add an IP address to a domain whitelabel.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

## URI Parameters
| URI Parameter   | Type  |  Description  |
|---|---|---|
| id | integer  | ID of the domain to which you are adding an IP |

### POST /whitelabel/domains/{id}/ips


```php
$request_body = json_decode('{
  "ip": "192.168.0.1"
}');
$id = "test_url_param";
$response = $sg->client->whitelabel()->domains()->_($id)->ips()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Remove an IP from a domain whitelabel.

**This endpoint allows you to remove a domain's IP address from that domain's whitelabel.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

## URI Parameters
| URI Parameter   | Type  | Description  |
|---|---|---|
| id | integer  | ID of the domain whitelabel to delete the IP from. |
| ip | string | IP to remove from the domain whitelabel. |

### DELETE /whitelabel/domains/{id}/ips/{ip}


```php
$id = "test_url_param";
$ip = "test_url_param";
$response = $sg->client->whitelabel()->domains()->_($id)->ips()->_($ip)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Validate a domain whitelabel.

**This endpoint allows you to validate a domain whitelabel. If it fails, it will return an error message describing why the whitelabel could not be validated.**

A domain whitelabel allows you to remove the via or sent on behalf of message that your recipients see when they read your emails. Whitelabeling a domain allows you to replace sendgrid.net with your personal sending domain. You will be required to create a subdomain so that SendGrid can generate the DNS records which you must give to your host provider. If you choose to use Automated Security, SendGrid will provide you with 3 CNAME records. If you turn Automated Security off, you will be given 2 TXT records and 1 MX record.

For more information on whitelabeling, please see our [User Guide](https://sendgrid.com/docs/User_Guide/Settings/Whitelabel/index.html)

## URI Parameters
| URI Parameter   | Type   | Description  |
|---|---|---|
| id | integer  |ID of the domain whitelabel to validate. |

### POST /whitelabel/domains/{id}/validate


```php
$id = "test_url_param";
$response = $sg->client->whitelabel()->domains()->_($id)->validate()->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create an IP whitelabel

**This endpoint allows you to create an IP whitelabel.**

When creating an IP whitelable, you should use the same subdomain that you used when you created a domain whitelabel.

A IP whitelabel consists of a subdomain and domain that will be used to generate a reverse DNS record for a given IP. Once SendGrid has verified that the appropriate A record for the IP has been created, the appropriate reverse DNS record for the IP is generated.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/ips.html).

### POST /whitelabel/ips


```php
$request_body = json_decode('{
  "domain": "example.com",
  "ip": "192.168.1.1",
  "subdomain": "email"
}');
$response = $sg->client->whitelabel()->ips()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all IP whitelabels

**This endpoint allows you to retrieve all of the IP whitelabels that have been createdy by this account.**

You may include a search key by using the "ip" parameter. This enables you to perform a prefix search for a given IP segment (e.g. "192.").

A IP whitelabel consists of a subdomain and domain that will be used to generate a reverse DNS record for a given IP. Once SendGrid has verified that the appropriate A record for the IP has been created, the appropriate reverse DNS record for the IP is generated.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/ips.html).

### GET /whitelabel/ips


```php
$query_params = json_decode('{"ip": "test_string", "limit": 1, "offset": 1}');
$response = $sg->client->whitelabel()->ips()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve an IP whitelabel

**This endpoint allows you to retrieve an IP whitelabel.**

A IP whitelabel consists of a subdomain and domain that will be used to generate a reverse DNS record for a given IP. Once SendGrid has verified that the appropriate A record for the IP has been created, the appropriate reverse DNS record for the IP is generated.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/ips.html).

### GET /whitelabel/ips/{id}


```php
$id = "test_url_param";
$response = $sg->client->whitelabel()->ips()->_($id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete an IP whitelabel

**This endpoint allows you to delete an IP whitelabel.**

A IP whitelabel consists of a subdomain and domain that will be used to generate a reverse DNS record for a given IP. Once SendGrid has verified that the appropriate A record for the IP has been created, the appropriate reverse DNS record for the IP is generated.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/ips.html).

### DELETE /whitelabel/ips/{id}


```php
$id = "test_url_param";
$response = $sg->client->whitelabel()->ips()->_($id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Validate an IP whitelabel

**This endpoint allows you to validate an IP whitelabel.**

A IP whitelabel consists of a subdomain and domain that will be used to generate a reverse DNS record for a given IP. Once SendGrid has verified that the appropriate A record for the IP has been created, the appropriate reverse DNS record for the IP is generated.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/ips.html).

### POST /whitelabel/ips/{id}/validate


```php
$id = "test_url_param";
$response = $sg->client->whitelabel()->ips()->_($id)->validate()->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Create a Link Whitelabel

**This endpoint allows you to create a new link whitelabel.**

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### POST /whitelabel/links


```php
$request_body = json_decode('{
  "default": true,
  "domain": "example.com",
  "subdomain": "mail"
}');
$query_params = json_decode('{"limit": 1, "offset": 1}');
$response = $sg->client->whitelabel()->links()->post($request_body, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve all link whitelabels

**This endpoint allows you to retrieve all link whitelabels.**

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### GET /whitelabel/links


```php
$query_params = json_decode('{"limit": 1}');
$response = $sg->client->whitelabel()->links()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a Default Link Whitelabel

**This endpoint allows you to retrieve the default link whitelabel.**

Default link whitelabel is the actual link whitelabel to be used when sending messages. If there are multiple link whitelabels, the default is determined by the following order:
<ul>
  <li>Validated link whitelabels marked as "default"</li>
  <li>Legacy link whitelabels (migrated from the whitelabel wizard)</li>
  <li>Default SendGrid link whitelabel (i.e. 100.ct.sendgrid.net)</li>
</ul>

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### GET /whitelabel/links/default


```php
$query_params = json_decode('{"domain": "test_string"}');
$response = $sg->client->whitelabel()->links()->default()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve Associated Link Whitelabel

**This endpoint allows you to retrieve the associated link whitelabel for a subuser.**

Link whitelables can be associated with subusers from the parent account. This functionality allows
subusers to send mail using their parent's linke whitelabels. To associate a link whitelabel, the parent account
must first create a whitelabel and validate it. The parent may then associate that whitelabel with a subuser via the API or the Subuser Management page in the user interface.

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### GET /whitelabel/links/subuser


```php
$query_params = json_decode('{"username": "test_string"}');
$response = $sg->client->whitelabel()->links()->subuser()->get(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Disassociate a Link Whitelabel

**This endpoint allows you to disassociate a link whitelabel from a subuser.**

Link whitelables can be associated with subusers from the parent account. This functionality allows
subusers to send mail using their parent's linke whitelabels. To associate a link whitelabel, the parent account
must first create a whitelabel and validate it. The parent may then associate that whitelabel with a subuser via the API or the Subuser Management page in the user interface.

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### DELETE /whitelabel/links/subuser


```php
$query_params = json_decode('{"username": "test_string"}');
$response = $sg->client->whitelabel()->links()->subuser()->delete(null, $query_params);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Update a Link Whitelabel

**This endpoint allows you to update a specific link whitelabel. You can use this endpoint to change a link whitelabel's default status.**

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### PATCH /whitelabel/links/{id}


```php
$request_body = json_decode('{
  "default": true
}');
$id = "test_url_param";
$response = $sg->client->whitelabel()->links()->_($id)->patch($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Retrieve a Link Whitelabel

**This endpoint allows you to retrieve a specific link whitelabel.**

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### GET /whitelabel/links/{id}


```php
$id = "test_url_param";
$response = $sg->client->whitelabel()->links()->_($id)->get();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Delete a Link Whitelabel

**This endpoint allows you to delete a link whitelabel.**

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### DELETE /whitelabel/links/{id}


```php
$id = "test_url_param";
$response = $sg->client->whitelabel()->links()->_($id)->delete();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Validate a Link Whitelabel

**This endpoint allows you to validate a link whitelabel.**

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### POST /whitelabel/links/{id}/validate


```php
$id = "test_url_param";
$response = $sg->client->whitelabel()->links()->_($id)->validate()->post();
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```
## Associate a Link Whitelabel

**This endpoint allows you to associate a link whitelabel with a subuser account.**

Link whitelables can be associated with subusers from the parent account. This functionality allows
subusers to send mail using their parent's linke whitelabels. To associate a link whitelabel, the parent account
must first create a whitelabel and validate it. The parent may then associate that whitelabel with a subuser via the API or the Subuser Management page in the user interface.

Email link whitelabels allow all of the click-tracked links you send in your emails to include the URL of your domain instead of sendgrid.net.

For more information, please see our [User Guide](https://sendgrid.com/docs/API_Reference/Web_API_v3/Whitelabel/links.html).

### POST /whitelabel/links/{link_id}/subuser


```php
$request_body = json_decode('{
  "username": "jane@example.com"
}');
$link_id = "test_url_param";
$response = $sg->client->whitelabel()->links()->_($link_id)->subuser()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```

If you have a non-library SendGrid issue, please contact our [support team](https://support.sendgrid.com).

If you can't find a solution below, please open an [issue](https://github.com/sendgrid/sendgrid-php/issues).


## Table of Contents

* [Migrating from v2 to v3](#migrating)
* [Continue Using v2](#v2)
* [Testing v3 /mail/send Calls Directly](#testing)
* [Error Messages](#error)
* [Versions](#versions)
* [Environment Variables and Your SendGrid API Key](#environment)
* [Using the Package Manager](#package-manager)

<a name="migrating"></a>
## Migrating from v2 to v3

Please review [our guide](https://sendgrid.com/docs/Classroom/Send/v3_Mail_Send/how_to_migrate_from_v2_to_v3_mail_send.html) on how to migrate from v2 to v3.

<a name="v2"></a>
## Continue Using v2

[Here](https://github.com/sendgrid/sendgrid-php/tree/75970eb82f5629e66db4d6da08ff7ef0c507e9b0) is the last working version with v2 support.

Using composer:

```json
{
  "require": {
    "sendgrid/sendgrid": "~3.2"
  }
}
```

Download packaged zip [here](https://sendgrid-open-source.s3.amazonaws.com/sendgrid-php/versions/sendgrid-php-75970eb.zip).

<a name="testing"></a>
## Testing v3 /mail/send Calls Directly

[Here](https://sendgrid.com/docs/Classroom/Send/v3_Mail_Send/curl_examples.html) are some cURL examples for common use cases.

<a name="error"></a>
## Error Messages

To read the error message returned by SendGrid's API:

```php
try {
    $response = $sendgrid->client->mail()->send()->post($mail);
} catch (Exception $e) {
    echo 'Caught exception: ',  $e->getMessage(), "\n";
}
```

<a name="versions"></a>
## Versions

We follow the MAJOR.MINOR.PATCH versioning scheme as described by [SemVer.org](http://semver.org). Therefore, we recommend that you always pin (or vendor) the particular version you are working with to your code and never auto-update to the latest version. Especially when there is a MAJOR point release, since that is guarenteed to be a breaking change. Changes are documented in the [CHANGELOG](https://github.com/sendgrid/sendgrid-php/blob/master/CHANGELOG.md) and [releases](https://github.com/sendgrid/sendgrid-php/releases) section.

<a name="environment"></a>
## Environment Variables and Your SendGrid API Key

All of our examples assume you are using [environment variables](https://github.com/sendgrid/sendgrid-php#setup-environment-variables) to hold your SendGrid API key.

If you choose to add your SendGrid API key directly (not recommended):

`$apiKey = getenv('SENDGRID_API_KEY');`

becomes

`$apiKey = 'SENDGRID_API_KEY';`

In the first case SENDGRID_API_KEY is in reference to the name of the environment variable, while the second case references the actual SendGrid API Key.

<a name="package-manager"></a>
## Using the Package Manager

We upload this library to [Packagist](https://packagist.org/packages/sendgrid/sendgrid) whenever we make a release. This allows you to use [composer](https://getcomposer.org) for easy installation.

In most cases we recommend you download the latest version of the library, but if you need a different version, please use:

```json
{
  "require": {
    "sendgrid/sendgrid": "~X.X.X"
  }
}
```[![BuildStatus](https://travis-ci.org/sendgrid/sendgrid-php.svg?branch=master)](https://travis-ci.org/sendgrid/sendgrid-php)

Please see our announcement regarding [breaking changes](https://github.com/sendgrid/sendgrid-php/issues/290). Your support is appreciated!

**This library allows you to quickly and easily use the SendGrid Web API v3 via PHP.**

Version 5.X.X of this library provides full support for all SendGrid [Web API v3](https://sendgrid.com/docs/API_Reference/Web_API_v3/index.html) endpoints, including the new [v3 /mail/send](https://sendgrid.com/blog/introducing-v3mailsend-sendgrids-new-mail-endpoint).

This library represents the beginning of a new path for SendGrid. We want this library to be community driven and SendGrid led. We need your help to realize this goal. To help make sure we are building the right things in the right order, we ask that you create [issues](https://github.com/sendgrid/sendgrid-php/issues) and [pull requests](https://github.com/sendgrid/sendgrid-php/blob/master/CONTRIBUTING.md) or simply upvote or comment on existing issues or pull requests.

Please browse the rest of this README for further detail.

We appreciate your continued support, thank you!

# Table of Contents

* [Installation](#installation)
* [Quick Start](#quick_start)
* [Usage](#usage)
* [Use Cases](#use_cases)
* [Announcements](#announcements)
* [Roadmap](#roadmap)
* [How to Contribute](#contribute)
* [Troubleshooting](#troubleshooting)
* [About](#about)

<a name="installation"></a>
# Installation

## Prerequisites

- PHP version 5.6 or 7.0
- The SendGrid service, starting at the [free level](https://sendgrid.com/free?source=sendgrid-php)

## Setup Environment Variables

Update the development environment with your [SENDGRID_API_KEY](https://app.sendgrid.com/settings/api_keys), for example:

```bash
echo "export SENDGRID_API_KEY='YOUR_API_KEY'" > sendgrid.env
echo "sendgrid.env" >> .gitignore
source ./sendgrid.env
```

## Install Package

Add SendGrid to your `composer.json` file. If you are not using [Composer](http://getcomposer.org), you should be. It's an excellent way to manage dependencies in your PHP application.

```json
{
  "require": {
    "sendgrid/sendgrid": "~5.0.9"
  }
}
```

Then at the top of your PHP script require the autoloader:

```bash
require 'vendor/autoload.php';
```

#### Alternative: Install package from zip

If you are not using Composer, simply download and install the **[latest packaged release of the library as a zip](https://sendgrid-open-source.s3.amazonaws.com/sendgrid-php/sendgrid-php.zip)**.

[**⬇︎ Download Packaged Library ⬇︎**](https://sendgrid-open-source.s3.amazonaws.com/sendgrid-php/sendgrid-php.zip)

Then require the library from package:

```php
require("path/to/sendgrid-php/sendgrid-php.php");
```

Previous versions of the library can be found in the [version index](https://sendgrid-open-source.s3.amazonaws.com/index.html).

## Dependencies

- The SendGrid Service, starting at the [free level](https://sendgrid.com/free?source=sendgrid-php)
- [php-HTTP-Client](https://github.com/sendgrid/php-http-client)

<a name="quick_start"></a>
# Quick Start

## Hello Email

The following is the minimum needed code to send an email with the [/mail/send Helper](https://github.com/sendgrid/sendgrid-php/tree/master/lib/helpers/mail) ([here](https://github.com/sendgrid/sendgrid-php/blob/master/examples/helpers/mail/example.php#L22) is a full example):

```php
<?php
// If you are using Composer
require 'vendor/autoload.php';

// If you are not using Composer (recommended)
// require("path/to/sendgrid-php/sendgrid-php.php");

$from = new SendGrid\Email(null, "test@example.com");
$subject = "Hello World from the SendGrid PHP Library!";
$to = new SendGrid\Email(null, "test@example.com");
$content = new SendGrid\Content("text/plain", "Hello, Email!");
$mail = new SendGrid\Mail($from, $subject, $to, $content);

$apiKey = getenv('SENDGRID_API_KEY');
$sg = new \SendGrid($apiKey);

$response = $sg->client->mail()->send()->post($mail);
echo $response->statusCode();
echo $response->headers();
echo $response->body();
```

The `SendGrid\Mail` constructor creates a [personalization object](https://sendgrid.com/docs/Classroom/Send/v3_Mail_Send/personalizations.html) for you. [Here](https://github.com/sendgrid/sendgrid-php/blob/master/examples/helpers/mail/example.php#L16) is an example of how to add to it.

### Without Mail Helper Class

The following is the minimum needed code to send an email without the /mail/send Helper ([here](https://github.com/sendgrid/sendgrid-php/blob/master/examples/mail/mail.php#L28) is a full example):

```php
<?php
// If you are using Composer
require 'vendor/autoload.php';

// If you are not using Composer (recommended)
// require("path/to/sendgrid-php/sendgrid-php.php");

$request_body = json_decode('{
  "personalizations": [
    {
      "to": [
        {
          "email": "test@example.com"
        }
      ],
      "subject": "Hello World from the SendGrid PHP Library!"
    }
  ],
  "from": {
    "email": "test@example.com"
  },
  "content": [
    {
      "type": "text/plain",
      "value": "Hello, Email!"
    }
  ]
}');

$apiKey = getenv('SENDGRID_API_KEY');
$sg = new \SendGrid($apiKey);

$response = $sg->client->mail()->send()->post($request_body);
echo $response->statusCode();
echo $response->body();
echo $response->headers();
```

## General v3 Web API Usage (With Fluent Interface)

```php
<?php
// If you are using Composer (recommended)
require 'vendor/autoload.php';

// If you are not using Composer
// require("path/to/sendgrid-php/sendgrid-php.php");

$apiKey = getenv('SENDGRID_API_KEY');
$sg = new \SendGrid($apiKey);

$response = $sg->client->suppressions()->bounces()->get();

print $response->statusCode();
print $response->headers();
print $response->body();
```

## General v3 Web API Usage (Without Fluent Interface)

```php
<?php
// If you are using Composer (recommended)
require 'vendor/autoload.php';

// If you are not using Composer
// require("path/to/sendgrid-php/sendgrid-php.php");

$apiKey = getenv('SENDGRID_API_KEY');
$sg = new \SendGrid($apiKey);

$response = $sg->client->_("suppression/bounces")->get();

print $response->statusCode();
print $response->headers();
print $response->body();
```

<a name="usage"></a>
# Usage

- [SendGrid Docs](https://sendgrid.com/docs/API_Reference/index.html)
- [Library Usage
    Documentation](https://github.com/sendgrid/sendgrid-php/tree/master/USAGE.md)
- [Example Code](https://github.com/sendgrid/sendgrid-php/tree/master/examples)
- [How-to: Migration from v2 to v3](https://sendgrid.com/docs/Classroom/Send/v3_Mail_Send/how_to_migrate_from_v2_to_v3_mail_send.html)
- [v3 Web API Mail Send Helper](https://github.com/sendgrid/sendgrid-php/tree/master/lib/helpers/mail/README.md) - build a request object payload for a v3 /mail/send API call.

<a name="use_cases">
# Use Cases

[Examples of common API use cases](https://github.com/sendgrid/sendgrid-php/blob/master/USE_CASES.md), such as how to send an email with a transactional template.

<a name="announcements"></a>
# Announcements

Please see our announcement regarding [breaking changes](https://github.com/sendgrid/sendgrid-php/issues/290). Your support is appreciated!

All updates to this library is documented in our [CHANGELOG](https://github.com/sendgrid/sendgrid-php/blob/master/CHANGELOG.md) and [releases](https://github.com/sendgrid/sendgrid-php/releases)

<a name="roadmap"></a>
# Roadmap

If you are interested in the future direction of this project, please take a look at our open [issues](https://github.com/sendgrid/sendgrid-php/issues) and [pull requests](https://github.com/sendgrid/sendgrid-php/pulls). We would love to hear your feedback.

<a name="contribute"></a>
# How to Contribute

We encourage contribution to our libraries (you might even score some nifty swag), please see our [CONTRIBUTING](https://github.com/sendgrid/sendgrid-php/blob/master/CONTRIBUTING.md) guide for details.

Quick links:

- [Feature Request](https://github.com/sendgrid/sendgrid-php/blob/master/CONTRIBUTING.md#feature_request)
- [Bug Reports](https://github.com/sendgrid/sendgrid-php/blob/master/CONTRIBUTING.md#submit_a_bug_report)
- [Sign the CLA to Create a Pull Request](https://github.com/sendgrid/sendgrid-php/blob/master/CONTRIBUTING.md#cla)
- [Improvements to the Codebase](https://github.com/sendgrid/sendgrid-php/blob/master/CONTRIBUTING.md#improvements_to_the_codebase)

<a name="troubleshooting"></a>
# Troubleshooting

Please see our [troubleshooting guide](https://github.com/sendgrid/sendgrid-php/blob/master/TROUBLESHOOTING.md) for common library issues.

<a name="about"></a>
# About

sendgrid-php is guided and supported by the SendGrid [Developer Experience Team](mailto:dx@sendgrid.com).

sendgrid-php is maintained and funded by SendGrid, Inc. The names and logos for sendgrid-php are trademarks of SendGrid, Inc.

![SendGrid Logo]
(https://uiux.s3.amazonaws.com/2016-logos/email-logo%402x.png)

This documentation provides examples for specific use cases. Please [open an issue](https://github.com/sendgrid/sendgrid-php/issues) or make a pull request for any use cases you would like us to document here. Thank you!

# Table of Contents

* [Transactional Templates](#transactional_templates)

<a name="transactional_templates"></a>
# Transactional Templates

For this example, we assume you have created a [transactional template](https://sendgrid.com/docs/User_Guide/Transactional_Templates/index.html). Following is the template content we used for testing.

Template ID (replace with your own):

```text
13b8f94f-bcae-4ec6-b752-70d6cb59f932
```

Email Subject:

```text
<%subject%>
```

Template Body:

```html
<html>
<head>
	<title></title>
</head>
<body>
Hello -name-,
<br /><br/>
I'm glad you are trying out the template feature!
<br /><br/>
<%body%>
<br /><br/>
I hope you are having a great day in -city- :)
<br /><br/>
</body>
</html>
```

## With Mail Helper Class

```php
<?php
// If you are using Composer
require 'vendor/autoload.php';

// If you are not using Composer (recommended)
// require("path/to/sendgrid-php/sendgrid-php.php");

$from = new SendGrid\Email(null, "test@example.com");
$subject = "I'm replacing the subject tag";
$to = new SendGrid\Email(null, "test@example.com");
$content = new SendGrid\Content("text/html", "I'm replacing the <strong>body tag</strong>");
$mail = new SendGrid\Mail($from, $subject, $to, $content);
$mail->personalization[0]->addSubstitution("-name-", "Example User");
$mail->personalization[0]->addSubstitution("-city-", "Denver");
$mail->setTemplateId("13b8f94f-bcae-4ec6-b752-70d6cb59f932");

$apiKey = getenv('SENDGRID_API_KEY');
$sg = new \SendGrid($apiKey);

try {
    $response = $sg->client->mail()->send()->post($mail);
} catch (Exception $e) {
    echo 'Caught exception: ',  $e->getMessage(), "\n";
}

echo $response->statusCode();
echo $response->headers();
echo $response->body();
```

## Without Mail Helper Class

```php
<?php
// If you are using Composer
require 'vendor/autoload.php';

// If you are not using Composer (recommended)
// require("path/to/sendgrid-php/sendgrid-php.php");

$request_body = json_decode('{
  "personalizations": [
    {
      "to": [
        {
          "email": "dx@sendgrid.com"
        }
      ],
      "substitutions": {
        "-name-": "Example User",
        "-city-": "Denver"
      },
      "subject": "I\'m replacing the subject tag"
    }
  ],
  "from": {
    "email": "elmer@sendgrid.com"
  },
  "content": [
    {
      "type": "text/html",
      "value": "I\'m replacing the <strong>body tag</strong>"
    }
  ],
  "template_id": "13b8f94f-bcae-4ec6-b752-70d6cb59f932"
}');

$apiKey = getenv('SENDGRID_API_KEY');
$sg = new \SendGrid($apiKey);

try {
    $response = $sg->client->mail()->send()->post($request_body);
} catch (Exception $e) {
    echo 'Caught exception: ',  $e->getMessage(), "\n";
}

echo $response->statusCode();
echo $response->body();
echo $response->headers();
```# Change Log
All notable changes to this project will be documented in this file.

This project adheres to [Semantic Versioning](http://semver.org/).

## [5.0.9] - 2016-09-13 ##
### Fixed
- Pull request #289: [Replace "\jsonSerializable" with "\JsonSerializable" ](https://github.com/sendgrid/sendgrid-php/pull/289)
- Thanks to [Issei.M](https://github.com/issei-m) for the pull request!

## [5.0.8] - 2016-08-24 ##
### Added
- Table of Contents in the README
- Added a [USE_CASES.md](https://github.com/sendgrid/sendgrid-php/blob/master/USE_CASES.md) section, with the first use case example for transactional templates

## [5.0.7] - 2016-07-25 ##
### Added
- [Troubleshooting](https://github.com/sendgrid/sendgrid-php/blob/master/TROUBLESHOOTING.md) section

## [5.0.6] - 2016-07-20 ##
### Added
- README updates
- Update introduction blurb to include information regarding our forward path
- Update the v3 /mail/send example to include non-helper usage
- Update the generic v3 example to include non-fluent interface usage

## [5.0.5] - 2016-07-12 ##
### Added
- Update docs, unit tests and examples to include Sender ID

## [5.0.4] - 2016-07-07 ##
### Added
- Tests now mocked automatically against [prism](https://stoplight.io/prism/)

## [5.0.3] - 2016-07-05 ##
### Updated
- Content based on our updated [Swagger/OAI doc](https://github.com/sendgrid/sendgrid-oai)

## [5.0.2] - 2016-07-05 ##
### Added
- Accept: application/json header per https://sendgrid.com/docs/API_Reference/Web_API_v3/How_To_Use_The_Web_API_v3/requests.html

### Updated
- Content based on our updated [Swagger/OAI doc](https://github.com/sendgrid/sendgrid-oai)

## [5.0.1] - 2016-06-17 ##
### Fixed
- Issue with packaged version for non-composer uses

## [5.0.0] - 2016-06-13 ##
### Added
- Breaking change to support the v3 Web API
- New HTTP client
- v3 Mail Send helper

## [v4.0.4] - (2016-02-18) ##
### Added
- Ability to add scopes to API Keys endpoint [POST]

## [v4.0.3] - (2016-02-18) ##
### Added
- API Keys endpoint [PUT]

## [v4.0.2] - (2015-12-15) ##
### Added
- Tests for API Keys endpoint [POST, PATCH, DELETE]

## [v4.0.1] - (2015-12-03) ##
### Fixed
- HTTP 406 Not Acceptable Errors [#177](https://github.com/sendgrid/sendgrid-php/issues/177)

## [v4.0.0] - (2015-10-16) ##
### Added
- Added support for accessing the [SendGrid Web API v3 endpoints](https://sendgrid.com/docs/API_Reference/Web_API_v3/index.html)
- Implemented part of the /api_keys, /groups and /suppressions endpoints

## [v3.2.0] - (2015-05-13) ##

### Added
- Specify Guzzle proxy via [#149](https://github.com/sendgrid/sendgrid-php/pull/149)
- Option to disable exception raising

## [v3.1.0] - (2015-04-27)
### Added
- Support for API keys

## [v3.0.0] - (2015-04-14)
### Fixed
- CC and BCC not working with SMTPAPI To

### Changed
- **Breaking:** A `\SendGrid\Exception` is now raised when response is not 200
- **Breaking:** `addTo` now uses the Web API parameter as opposed to the SMTPAPI Header. Substitutions will most likely break unless you update to use `addSmtpapiTo`
- The library now depends on Guzzle3
- Major refactoring

### Added
- **Breaking:** `send()` now returns an instance of `\SendGrid\Response`
- Numerous missing methods for new functionality
- `addSmtpapiTo` for using the SMTPAPI To

## [v2.2.1] - (2014-01-29)
### Fixed
- Fix turn_off_ssl_verification option via [#123](https://github.com/sendgrid/sendgrid-php/pull/123)

## [v2.2.0] - (2014-01-12)
### Changed
- Remove [Unirest](https://github.com/Mashape/unirest-php/) and replace with native cURL
- Add CHANGELOG.md
Hello! Thank you for choosing to help contribute to one of the SendGrid open source libraries. There are many ways you can contribute and help is always welcome.  We simply ask that you follow the following contribution policies.

- [CLAs and CCLAs](#cla)
- [Roadmap & Milestones](#roadmap)
- [Feature Request](#feature_request)
- [Submit a Bug Report](#submit_a_bug_report)
- [Improvements to the Codebase](#improvements_to_the_codebase)
- [Understanding the Code Base](#understanding_the_codebase)
- [Testing](#testing)
- [Style Guidelines & Naming Conventions](#style_guidelines_and_naming_conventions)
- [Creating a Pull Request](#creating_a_pull_request)

<a name="roadmap"></a>
We use [Milestones](https://github.com/sendgrid/sendgrid-php/milestones) to help define current roadmaps, please feel free to grab an issue from the current milestone. Please indicate that you have begun work on it to avoid collisions. Once a PR is made, community review, comments, suggestions and additional PRs are welcomed and encouraged.

<a name="cla"></a>
## CLAs and CCLAs

Before you get started, SendGrid requires that a SendGrid Contributor License Agreement (CLA) or a SendGrid Company Contributor Licensing Agreement (CCLA) be filled out by every contributor to a SendGrid open source project.

Our goal with the CLA and CCLA is to clarify the rights of our contributors and reduce other risks arising from inappropriate contributions.  The CLA also clarifies the rights SendGrid holds in each contribution and helps to avoid misunderstandings over what rights each contributor is required to grant to SendGrid when making a contribution.  In this way the CLA and CCLA encourage broad participation by our open source community and help us build strong open source projects, free from any individual contributor withholding or revoking rights to any contribution.

SendGrid does not merge a pull request made against a SendGrid open source project until that pull request is associated with a signed CLA (or CCLA). Copies of the CLA and CCLA are available [here](https://drive.google.com/a/sendgrid.com/file/d/0B0PlcM9qA91LN2VEUTJWU2RIVXc/view).

You may submit your completed [CLA or CCLA](https://drive.google.com/a/sendgrid.com/file/d/0B0PlcM9qA91LN2VEUTJWU2RIVXc/view) to SendGrid at [dx@sendgrid.com](mailto:dx@sendgrid.com).  SendGrid will then confirm you are ready to begin making contributions.

There are a few ways to contribute, which we'll enumerate below:

<a name="feature_request"></a>
## Feature Request

If you'd like to make a feature request, please read this section.

The GitHub issue tracker is the preferred channel for library feature requests, but please respect the following restrictions:

- Please **search for existing issues** in order to ensure we don't have duplicate bugs/feature requests.
- Please be respectful and considerate of others when commenting on issues

<a name="submit_a_bug_report"></a>
## Submit a Bug Report

Note: DO NOT include your credentials in ANY code examples, descriptions, or media you make public.

A software bug is a demonstrable issue in the code base. In order for us to diagnose the issue and respond as quickly as possible, please add as much detail as possible into your bug report.

Before you decide to create a new issue, please try the following:

1. Check the Github issues tab if the identified issue has already been reported, if so, please add a +1 to the existing post.
2. Update to the latest version of this code and check if issue has already been fixed
3. Copy and fill in the Bug Report Template we have provided below

### Please use our Bug Report Template

In order to make the process easier, we've included a [sample bug report template](https://github.com/sendgrid/sendgrid-php/.github/ISSUE_TEMPLATE) (borrowed from [Ghost](https://github.com/TryGhost/Ghost/)). The template uses [GitHub flavored markdown](https://help.github.com/articles/github-flavored-markdown/) for formatting.

<a name="improvements_to_the_codebase"></a>
## Improvements to the Codebase

We welcome direct contributions to the sendgrid-php code base. Thank you!

### Development Environment ###

#### Install and Run Locally ####

##### Prerequisites #####

- PHP 5.6 or 7.0

##### Initial setup: #####

```bash
git clone https://github.com/sendgrid/sendgrid-php.git
cd sendgrid-php
composer install
```

## Environment Variables

First, get your free SendGrid account [here](https://sendgrid.com/free?source=sendgrid-php).

Next, update your environment with your [SENDGRID_API_KEY](https://app.sendgrid.com/settings/api_keys).

```bash
echo "export SENDGRID_API_KEY='YOUR_API_KEY'" > sendgrid.env
echo "sendgrid.env" >> .gitignore
source ./sendgrid.env
```

##### Execute: #####

See the [examples folder](https://github.com/sendgrid/sendgrid-php/tree/master/examples) to get started quickly.

If you are using composer, replace <PATH_TO> with the path to your `vendor/autoload.php`. Otherwise, include lib/SendGrid.php and lib/helpers/mail/Mail.php.

<a name="understanding_the_codebase"></a>
## Understanding the Code Base

**/examples**

Working examples that demonstrate usage.

```bash
php examples/example.php
```

**/test/unit**

Unit tests for the HTTP client.

**/lib**

The interface to the SendGrid API.

<a name="testing"></a>
## Testing

All PRs require passing tests before the PR will be reviewed.

All test files are in the [`/test/unit`](https://github.com/sendgrid/sendgrid-php/tree/master/test/unit) directory.

For the purposes of contributing to this repo, please update the [`SendGridTest.php`](https://github.com/sendgrid/sendgrid-php/tree/master/test/unit/SendGridTest.php) file with unit tests as you modify the code.

```bash
composer install
cd test/unit
../../vendor/bin/phpunit . --bootstrap bootstrap.php --filter test*
```

<a name="style_guidelines_and_naming_conventions"></a>
## Style Guidelines & Naming Conventions

Generally, we follow the style guidelines as suggested by the official language. However, we ask that you conform to the styles that already exist in the library. If you wish to deviate, please explain your reasoning.

- [pear coding standards](https://pear.php.net/manual/en/standards.php)

Please run your code through:

- [PHP Code Sniffer](https://github.com/squizlabs/PHP_CodeSniffer)

## Creating a Pull Request<a name="creating_a_pull_request"></a>

1. [Fork](https://help.github.com/fork-a-repo/) the project, clone your fork,
   and configure the remotes:

   ```bash
   # Clone your fork of the repo into the current directory
   git clone https://github.com/sendgrid/sendgrid-php
   # Navigate to the newly cloned directory
   cd sendgrid-python
   # Assign the original repo to a remote called "upstream"
   git remote add upstream https://github.com/sendgrid/sendgrid-php
   ```

2. If you cloned a while ago, get the latest changes from upstream:

   ```bash
   git checkout <dev-branch>
   git pull upstream <dev-branch>
   ```

3. Create a new topic branch (off the main project development branch) to
   contain your feature, change, or fix:

   ```bash
   git checkout -b <topic-branch-name>
   ```

4. Commit your changes in logical chunks. Please adhere to these [git commit
   message guidelines](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
   or your code is unlikely be merged into the main project. Use Git's
   [interactive rebase](https://help.github.com/articles/interactive-rebase)
   feature to tidy up your commits before making them public.

4a. Create tests.

4b. Create or update the example code that demonstrates the functionality of this change to the code.

5. Locally merge (or rebase) the upstream development branch into your topic branch:

   ```bash
   git pull [--rebase] upstream master
   ```

6. Push your topic branch up to your fork:

   ```bash
   git push origin <topic-branch-name>
   ```

7. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/)
    with a clear title and description against the `master` branch. All tests must be passing before we will review the PR.

If you have any additional questions, please feel free to [email](mailto:dx@sendgrid.com) us or create an issue in this repo.**This helper allows you to quickly and easily build a Mail object for sending email through SendGrid.**

# Quick Start

Run the [example](https://github.com/sendgrid/sendgrid-php/blob/master/examples/helpers/mail/example.php) (make sure you have set your environment variable to include your SENDGRID_API_KEY).

```bash
php examples/helpers/mail/example.php
```

## Usage

- See this complete working [example](https://github.com/sendgrid/sendgrid-php/blob/master/examples/helpers/mail/example.php).
- [Documentation](https://sendgrid.com/docs/API_Reference/Web_API_v3/Mail/overview.html)
[![Travis Badge](https://travis-ci.org/sendgrid/php-http-client.svg?branch=master)](https://travis-ci.org/sendgrid/php-http-client)

**Quickly and easily access any RESTful or RESTful-like API.**

If you are looking for the SendGrid API client library, please see [this repo](https://github.com/sendgrid/sendgrid-php).

# Announcements

All updates to this library is documented in our [CHANGELOG](https://github.com/sendgrid/php-http-client/blob/master/CHANGELOG.md).

# Installation

Add php-http-client to your `composer.json` file. If you are not using [Composer](http://getcomposer.org), you should be. It's an excellent way to manage dependencies in your PHP application.

```json
{
  "require": {
    "sendgrid/php-http-client": "3.*"
  }
}
```

Then at the top of your PHP script require the autoloader:

```php
require __DIR__ . '/vendor/autoload.php';
```

Then from the command line:

```bash
composer install
```

# Quick Start

Here is a quick example:

`GET /your/api/{param}/call`

```php
require 'vendor/autoload.php';
$global_headers = array(Authorization: Basic XXXXXXX);
$client = SendGrid\Client('base_url', 'global_headers');
$response = $client->your()->api()->_($param)->call()->get();
print $response->statusCode();
print $response->headers();
print $response->body();
```

`POST /your/api/{param}/call` with headers, query parameters and a request body with versioning.

```php
require 'vendor/autoload.php';
$global_headers = array(Authorization: Basic XXXXXXX);
$client = SendGrid\Client('base_url', 'global_headers');
$query_params = array('hello' => 0, 'world' => 1);
$request_headers = array('X-Test' => 'test');
$data = array('some' => 1, 'awesome' => 2, 'data' => 3);
$response = $client->your()->api()->_($param)->call()->post('data',
                                                            'query_params',
                                                            'request_headers');
print $response->statusCode();
print $response->headers();
print $response->body();
```

# Usage

- [Example Code](https://github.com/sendgrid/php-http-client/tree/master/examples)

## Roadmap

If you are intersted in the future direction of this project, please take a look at our [milestones](https://github.com/sendgrid/php-http-client/milestones). We would love to hear your feedback.

## How to Contribute

We encourage contribution to our libraries, please see our [CONTRIBUTING](https://github.com/sendgrid/php-http-client/blob/master/CONTRIBUTING.md)) guide for details.

Quick links:

- [Feature Request](https://github.com/sendgrid/php-http-client/blob/master/CONTRIBUTING.md#feature_request)
- [Bug Reports](https://github.com/sendgrid/php-http-client/blob/master/CONTRIBUTING.md#submit_a_bug_report)
- [Sign the CLA to Create a Pull Request](https://github.com/sendgrid/php-http-client/blob/master/CONTRIBUTING.md#cla)
- [Improvements to the Codebase](https://github.com/sendgrid/php-http-client/blob/master/CONTRIBUTING.md#improvements_to_the_codebase)

# Thanks

We were inspired by the work done on [birdy](https://github.com/inueni/birdy) and [universalclient](https://github.com/dgreisen/universalclient).

# About

php-http-client is guided and supported by the SendGrid [Developer Experience Team](mailto:dx@sendgrid.com).

php-http-client is maintained and funded by SendGrid, Inc. The names and logos for php-http-client are trademarks of SendGrid, Inc.

![SendGrid Logo]
(https://uiux.s3.amazonaws.com/2016-logos/email-logo%402x.png)
# Change Log
All notable changes to this project will be documented in this file.

This project adheres to [Semantic Versioning](http://semver.org/).

## [3.1.0] - 2016-06-10
### Added
- Automatically add Content-Type: application/json when there is a request body

## [3.0.0] - 2016-06-06
### Changed
- Made the Request and Response variables non-redundant. e.g. request.requestBody becomes request.body

## [2.0.2] - 2016-02-29
### Fixed
- Renaming files to conform to PSR-0, git ignored the case in 2.0.1

## [2.0.1] - 2016-02-29
### Fixed
- Renaming files to conform to PSR-0

## [1.0.1] - 2016-02-29
### Fixed
- Composer/Packagist install issues resolved

## [1.0.0] - 2016-02-29
### Added
- We are live!
Hello! Thank you for choosing to help contribute to one of the SendGrid open source projects. There are many ways you can contribute and help is always welcome.  We simply ask that you follow the following contribution policies.

- [CLAs and CCLAs](#cla)
- [Roadmap & Milestones](#roadmap)
- [Feature Request](#feature_request)
- [Submit a Bug Report](#submit_a_bug_report)
- [Improvements to the Codebase](#improvements_to_the_codebase)
- [Understanding the Code Base](#understanding_the_codebase)
- [Testing](#testing)
- [Style Guidelines & Naming Conventions](#style_guidelines_and_naming_conventions)
- [Creating a Pull Request](#creating_a_pull_request)

<a name="roadmap"></a>
We use [Milestones](https://github.com/sendgrid/php-http-client/milestones) to help define current roadmaps, please feel free to grab an issue from the current milestone. Please indicate that you have begun work on it to avoid collisions. Once a PR is made, community review, comments, suggestions and additional PRs are welcomed and encouraged.

<a name="cla"></a>
## CLAs and CCLAs

Before you get started, SendGrid requires that a SendGrid Contributor License Agreement (CLA) or a SendGrid Company Contributor Licensing Agreement (CCLA) be filled out by every contributor to a SendGrid open source project.

Our goal with the CLA and CCLA is to clarify the rights of our contributors and reduce other risks arising from inappropriate contributions.  The CLA also clarifies the rights SendGrid holds in each contribution and helps to avoid misunderstandings over what rights each contributor is required to grant to SendGrid when making a contribution.  In this way the CLA and CCLA encourage broad participation by our open source community and help us build strong open source projects, free from any individual contributor withholding or revoking rights to any contribution.

SendGrid does not merge a pull request made against a SendGrid open source project until that pull request is associated with a signed CLA (or CCLA). Copies of the CLA and CCLA are available [here](https://drive.google.com/a/sendgrid.com/file/d/0B0PlcM9qA91LN2VEUTJWU2RIVXc/view).

You may submit your completed [CLA or CCLA](https://drive.google.com/a/sendgrid.com/file/d/0B0PlcM9qA91LN2VEUTJWU2RIVXc/view) to SendGrid at [dx@sendgrid.com](mailto:dx@sendgrid.com).  SendGrid will then confirm you are ready to begin making contributions.

There are a few ways to contribute, which we'll enumerate below:

<a name="feature_request"></a>
## Feature Request

If you'd like to make a feature request, please read this section.

The GitHub issue tracker is the preferred channel for library feature requests, but please respect the following restrictions:

- Please **search for existing issues** in order to ensure we don't have duplicate bugs/feature requests.
- Please be respectful and considerate of others when commenting on issues

<a name="submit_a_bug_report"></a>
## Submit a Bug Report

Note: DO NOT include your credentials in ANY code examples, descriptions, or media you make public.

A software bug is a demonstrable issue in the code base. In order for us to diagnose the issue and respond as quickly as possible, please add as much detail as possible into your bug report.

Before you decide to create a new issue, please try the following:

1. Check the Github issues tab if the identified issue has already been reported, if so, please add a +1 to the existing post.
2. Update to the latest version of this code and check if issue has already been fixed
3. Copy and fill in the Bug Report Template we have provided below

### Please use our Bug Report Template

In order to make the process easier, we've included a [sample bug report template](https://github.com/sendgrid/php-http-client/.github/ISSUE_TEMPLATE) (borrowed from [Ghost](https://github.com/TryGhost/Ghost/)). The template uses [GitHub flavored markdown](https://help.github.com/articles/github-flavored-markdown/) for formatting.

<a name="improvements_to_the_codebase"></a>
## Improvements to the Codebase

We welcome direct contributions to the php-http-client code base. Thank you!

### Development Environment ###

#### Install and Run Locally ####

##### Prerequisites #####

- PHP 5.2 through 5.6
- [Composer](https://getcomposer.org/)

##### Initial setup: #####

```bash
git clone https://github.com/sendgrid/php-http-client.git
cd php-http-client
```

## Environment Variables

First, get your free SendGrid account [here](https://sendgrid.com/free?source=php-http-client).

Next, update your environment with your [SENDGRID_API_KEY](https://app.sendgrid.com/settings/api_keys).

```bash
echo "export SENDGRID_API_KEY='YOUR_API_KEY'" > sendgrid.env
echo "sendgrid.env" >> .gitignore
source ./sendgrid.env
```

##### Execute: #####

See the [examples folder](https://github.com/sendgrid/php-http-client/tree/master/examples
<a name="understanding_the_codebase"></a>
## Understanding the Code Base

**/examples**

Working examples that demonstrate usage.

**/test/unit**

Unit tests.

**/lib/SendGrid/Client.php**

An HTTP client with a fluent interface using method chaining and reflection. By returning self on [__call](https://github.com/sendgrid/php-http-client/blob/master/lib/client.php#L212) and [_()](https://github.com/sendgrid/php-http-client/blob/master/lib/client.pph#L198), we can dynamically build the URL using method chaining and [__call](https://github.com/sendgrid/php-http-client/blob/master/lib/client.php#L212) allows us to dynamically receive the method calls to achieve reflection.

This allows for the following mapping from a URL to a method chain:

`/api_client/{api_key_id}/version` maps to `client->api_client().->_($api_key_id)->version-><method>()` where <method> is a [HTTP verb](https://github.com/sendgrid/php-http-client/blob/master/lib/client.php#L94).

**/lib/SendGrid/Config.php**

Loads the environment variables.

<a name="testing"></a>
## Testing

All PRs require passing tests before the PR will be reviewed.

All test files are in the [`/test/unit`](https://github.com/sendgrid/php-http-client/tree/master/test/unit) directory.

For the purposes of contributing to this repo, please update the [`ClientTest.php`](https://github.com/sendgrid/php-http-client/blob/master/test/unit/ClientTest.php) file with unit tests as you modify the code.

```bash
phpunit --bootstrap test/unit/bootstrap.php --filter test* test/unit
```

<a name="style_guidelines_and_naming_conventions"></a>
## Style Guidelines & Naming Conventions

Generally, we follow the style guidelines as suggested by the official language. However, we ask that you conform to the styles that already exist in the library. If you wish to deviate, please explain your reasoning.

- [pear coding standards](https://pear.php.net/manual/en/standards.php)

Please run your code through:

- [PHP Code Sniffer](https://github.com/squizlabs/PHP_CodeSniffer)

## Creating a Pull Request<a name="creating_a_pull_request"></a>

1. [Fork](https://help.github.com/fork-a-repo/) the project, clone your fork,
   and configure the remotes:

   ```bash
   # Clone your fork of the repo into the current directory
   git clone https://github.com/sendgrid/php-http-client
   # Navigate to the newly cloned directory
   cd sendgrid-python
   # Assign the original repo to a remote called "upstream"
   git remote add upstream https://github.com/sendgrid/php-http-client
   ```

2. If you cloned a while ago, get the latest changes from upstream:

   ```bash
   git checkout <dev-branch>
   git pull upstream <dev-branch>
   ```

3. Create a new topic branch (off the main project development branch) to
   contain your feature, change, or fix:

   ```bash
   git checkout -b <topic-branch-name>
   ```

4. Commit your changes in logical chunks. Please adhere to these [git commit
   message guidelines](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
   or your code is unlikely be merged into the main project. Use Git's
   [interactive rebase](https://help.github.com/articles/interactive-rebase)
   feature to tidy up your commits before making them public.

4a. Create tests.

4b. Create or update the example code that demonstrates the functionality of this change to the code.

5. Locally merge (or rebase) the upstream development branch into your topic branch:

   ```bash
   git pull [--rebase] upstream master
   ```

6. Push your topic branch up to your fork:

   ```bash
   git push origin <topic-branch-name>
   ```

7. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/)
    with a clear title and description against the `master` branch. All tests must be passing before we will review the PR.

If you have any additional questions, please feel free to [email](mailto:dx@sendgrid.com) us or create an issue in this repo.# Font Awesome 5.0.4

Thanks for downloading Font Awesome! We're so excited you're here.

Our documentation is available online. Just head here:

https://fontawesome.com
# @fortawesome/fontawesome-free-brands - SVG with JavaScript version

> "I came here to chew bubblegum and install Font Awesome 5 - and I'm all out of bubblegum"

[![npm](https://img.shields.io/npm/v/@fortawesome/fontawesome-free-brands.svg?style=flat-square)](https://www.npmjs.com/package/@fortawesome/fontawesome-free-brands)

## Installation

```
$ npm i --save @fortawesome/fontawesome-free-brands
```

Or

```
$ yarn add @fortawesome/fontawesome-free-brands
```

## Documentation

Get started [here](https://fontawesome.com/get-started/svg-with-js). Continue your journey [here](https://fontawesome.com/how-to-use/svg-with-js).

Or go straight to the [API documentation](https://fontawesome.com/how-to-use/font-awesome-api).

## Issues and support

Start with [GitHub issues](https://github.com/FortAwesome/Font-Awesome/issues) and ping us on [Twitter](https://twitter.com/fontawesome) if you need to.
# @fortawesome/fontawesome-free-webfonts - Web Fonts with CSS version

> "I came here to chew bubblegum and install Font Awesome 5 - and I'm all out of bubblegum"

[![npm](https://img.shields.io/npm/v/@fortawesome/fontawesome-free-webfonts.svg?style=flat-square)](https://www.npmjs.com/package/@fortawesome/fontawesome-free-webfonts)

## Installation

```
$ npm i --save @fortawesome/fontawesome-free-webfonts
```

Or

```
$ yarn add @fortawesome/fontawesome-free-webfonts
```

## What is this package?

This package includes CSS, Less, SCSS, and the individual web font files needed to get Font Awesome 5 Web Fonts with CSS working in your freeject.

Inside `node_modules/@fortawesome/fontawesome-free-webfonts` you'll find the following directories:

* `css`
* `less`
* `scss`
* `webfonts`

It's up to you to configure your build tool to use these.

## Documentation

Get started [here](https://fontawesome.com/get-started/web-fonts-with-css). Continue your journey [here](https://fontawesome.com/how-to-use/web-fonts-with-css).

## Issues and support

Start with [GitHub issues](https://github.com/FortAwesome/Font-Awesome/issues) and ping us on [Twitter](https://twitter.com/fontawesome) if you need to.
# @fortawesome/fontawesome-free-regular - SVG with JavaScript version

> "I came here to chew bubblegum and install Font Awesome 5 - and I'm all out of bubblegum"

[![npm](https://img.shields.io/npm/v/@fortawesome/fontawesome-free-regular.svg?style=flat-square)](https://www.npmjs.com/package/@fortawesome/fontawesome-free-regular)

## Installation

```
$ npm i --save @fortawesome/fontawesome-free-regular
```

Or

```
$ yarn add @fortawesome/fontawesome-free-regular
```

## Documentation

Get started [here](https://fontawesome.com/get-started/svg-with-js). Continue your journey [here](https://fontawesome.com/how-to-use/svg-with-js).

Or go straight to the [API documentation](https://fontawesome.com/how-to-use/font-awesome-api).

## Issues and support

Start with [GitHub issues](https://github.com/FortAwesome/Font-Awesome/issues) and ping us on [Twitter](https://twitter.com/fontawesome) if you need to.
# @fortawesome/fontawesome-free-solid - SVG with JavaScript version

> "I came here to chew bubblegum and install Font Awesome 5 - and I'm all out of bubblegum"

[![npm](https://img.shields.io/npm/v/@fortawesome/fontawesome-free-solid.svg?style=flat-square)](https://www.npmjs.com/package/@fortawesome/fontawesome-free-solid)

## Installation

```
$ npm i --save @fortawesome/fontawesome-free-solid
```

Or

```
$ yarn add @fortawesome/fontawesome-free-solid
```

## Documentation

Get started [here](https://fontawesome.com/get-started/svg-with-js). Continue your journey [here](https://fontawesome.com/how-to-use/svg-with-js).

Or go straight to the [API documentation](https://fontawesome.com/how-to-use/font-awesome-api).

## Issues and support

Start with [GitHub issues](https://github.com/FortAwesome/Font-Awesome/issues) and ping us on [Twitter](https://twitter.com/fontawesome) if you need to.
# @fortawesome/fontawesome-common-types - SVG with JavaScript

> "I came here to chew bubblegum and install Font Awesome 5 - and I'm all out of bubblegum"

[![npm](https://img.shields.io/npm/v/@fortawesome/fontawesome-common-types.svg?style=flat-square)](https://www.npmjs.com/package/@fortawesome/fontawesome-common-types)

## What is this package?

Font Awesome 5 JavaScript packages support TypeScript. This package abstracts out some of the common definitions that those packages use.

## Here be dragons

If you are trying to import types from this package we *highly* recommend you do the following instead as *all types in this package are re-exported to the main fontawesome package*.

your.ts

```
import {
  IconName
} from `@fortawesome/fontawesome`

const myIcon: IconName = "..."
```

## Issues and support

Start with [GitHub issues](https://github.com/FortAwesome/Font-Awesome/issues) and ping us on [Twitter](https://twitter.com/fontawesome) if you need to.
# @fortawesome/fontawesome - SVG with JavaScript version

> "I came here to chew bubblegum and install Font Awesome 5 - and I'm all out of bubblegum"

[![npm](https://img.shields.io/npm/v/@fortawesome/fontawesome.svg?style=flat-square)](https://www.npmjs.com/package/@fortawesome/fontawesome)

## Installation

```
$ npm i --save @fortawesome/fontawesome
```

Or

```
$ yarn add @fortawesome/fontawesome
```

## Documentation

Get started [here](https://fontawesome.com/get-started/svg-with-js). Continue your journey [here](https://fontawesome.com/how-to-use/svg-with-js).

Or go straight to the [API documentation](https://fontawesome.com/how-to-use/font-awesome-api).

## Issues and support

Start with [GitHub issues](https://github.com/FortAwesome/Font-Awesome/issues) and ping us on [Twitter](https://twitter.com/fontawesome) if you need to.
![PHPMailer](https://raw.github.com/PHPMailer/PHPMailer/master/examples/images/phpmailer.png)

# PHPMailer - A full-featured email creation and transfer class for PHP

Build status: [![Build Status](https://travis-ci.org/PHPMailer/PHPMailer.svg)](https://travis-ci.org/PHPMailer/PHPMailer)
[![Scrutinizer Quality Score](https://scrutinizer-ci.com/g/PHPMailer/PHPMailer/badges/quality-score.png?s=3758e21d279becdf847a557a56a3ed16dfec9d5d)](https://scrutinizer-ci.com/g/PHPMailer/PHPMailer/)
[![Code Coverage](https://scrutinizer-ci.com/g/PHPMailer/PHPMailer/badges/coverage.png?s=3fe6ca5fe8cd2cdf96285756e42932f7ca256962)](https://scrutinizer-ci.com/g/PHPMailer/PHPMailer/)

[![Latest Stable Version](https://poser.pugx.org/phpmailer/phpmailer/v/stable.svg)](https://packagist.org/packages/phpmailer/phpmailer) [![Total Downloads](https://poser.pugx.org/phpmailer/phpmailer/downloads)](https://packagist.org/packages/phpmailer/phpmailer) [![Latest Unstable Version](https://poser.pugx.org/phpmailer/phpmailer/v/unstable.svg)](https://packagist.org/packages/phpmailer/phpmailer) [![License](https://poser.pugx.org/phpmailer/phpmailer/license.svg)](https://packagist.org/packages/phpmailer/phpmailer)

## Class Features

- Probably the world's most popular code for sending email from PHP!
- Used by many open-source projects: WordPress, Drupal, 1CRM, SugarCRM, Yii, Joomla! and many more
- Integrated SMTP support - send without a local mail server
- Send emails with multiple TOs, CCs, BCCs and REPLY-TOs
- Multipart/alternative emails for mail clients that do not read HTML email
- Support for UTF-8 content and 8bit, base64, binary, and quoted-printable encodings
- SMTP authentication with LOGIN, PLAIN, NTLM, CRAM-MD5 and Google's XOAUTH2 mechanisms over SSL and TLS transports
- Error messages in 47 languages!
- DKIM and S/MIME signing support
- Compatible with PHP 5.0 and later
- Much more!

## Why you might need it

Many PHP developers utilize email in their code. The only PHP function that supports this is the `mail()` function. However, it does not provide any assistance for making use of popular features such as HTML-based emails and attachments.

Formatting email correctly is surprisingly difficult. There are myriad overlapping RFCs, requiring tight adherence to horribly complicated formatting and encoding rules - the vast majority of code that you'll find online that uses the `mail()` function directly is just plain wrong!
*Please* don't be tempted to do it yourself - if you don't use PHPMailer, there are many other excellent libraries that you should look at before rolling your own - try SwiftMailer, Zend_Mail, eZcomponents etc.

The PHP `mail()` function usually sends via a local mail server, typically fronted by a `sendmail` binary on Linux, BSD and OS X platforms, however, Windows usually doesn't include a local mail server; PHPMailer's integrated SMTP implementation allows email sending on Windows platforms without a local mail server.

## License

This software is distributed under the [LGPL 2.1](http://www.gnu.org/licenses/lgpl-2.1.html) license. Please read LICENSE for information on the
software availability and distribution.

## Installation & loading

PHPMailer is available via [Composer/Packagist](https://packagist.org/packages/phpmailer/phpmailer) (using semantic versioning), so just add this line to your `composer.json` file:

```json
"phpmailer/phpmailer": "~5.2"
```

or

```sh
composer require phpmailer/phpmailer
```

If you want to use the Gmail XOAUTH2 authentication class, you will also need to add a dependency on the `league/oauth2-client` package.

Alternatively, copy the contents of the PHPMailer folder into one of the `include_path` directories specified in your PHP configuration. If you don't speak git or just want a tarball, click the 'zip' button at the top of the page in GitHub.

If you're not using composer's autoloader, PHPMailer provides an SPL-compatible autoloader, and that is the preferred way of loading the library - just `require '/path/to/PHPMailerAutoload.php';` and everything should work. The autoloader does not throw errors if it can't find classes so it prepends itself to the SPL list, allowing your own (or your framework's) autoloader to catch errors. SPL autoloading was introduced in PHP 5.1.0, so if you are using a version older than that you will need to require/include each class manually.

PHPMailer does *not* declare a namespace because namespaces were only introduced in PHP 5.3.

If you want to use Google's XOAUTH2 authentication mechanism, you need to be running at least PHP 5.4, and load the dependencies listed in `composer.json`.

### Minimal installation

While installing the entire package manually or with composer is simple, convenient and reliable, you may want to include only vital files in your project. At the very least you will need [class.phpmailer.php](https://github.com/PHPMailer/PHPMailer/tree/master/class.phpmailer.php). If you're using SMTP, you'll need [class.smtp.php](https://github.com/PHPMailer/PHPMailer/tree/master/class.smtp.php), and if you're using POP-before SMTP, you'll need [class.pop3.php](https://github.com/PHPMailer/PHPMailer/tree/master/class.pop3.php). For all of these, we recommend you use [the autoloader](https://github.com/PHPMailer/PHPMailer/tree/master/PHPMailerAutoload.php) too as otherwise you will either have to `require` all classes manually or use some other autoloader. You can skip the [language](https://github.com/PHPMailer/PHPMailer/tree/master/language/) folder if you're not showing errors to users and can make do with English-only errors. You may need the additional classes in the [extras](extras/) folder if you are using those features, including NTLM authentication and ics generation. If you're using Google XOAUTH2 you will need `class.phpmaileroauth.php` and `class.oauth.php` classes too, as well as the composer dependencies.

## A Simple Example

```php
<?php
require 'PHPMailerAutoload.php';

$mail = new PHPMailer;

//$mail->SMTPDebug = 3;                               // Enable verbose debug output

$mail->isSMTP();                                      // Set mailer to use SMTP
$mail->Host = 'smtp1.example.com;smtp2.example.com';  // Specify main and backup SMTP servers
$mail->SMTPAuth = true;                               // Enable SMTP authentication
$mail->Username = 'user@example.com';                 // SMTP username
$mail->Password = 'secret';                           // SMTP password
$mail->SMTPSecure = 'tls';                            // Enable TLS encryption, `ssl` also accepted
$mail->Port = 587;                                    // TCP port to connect to

$mail->setFrom('from@example.com', 'Mailer');
$mail->addAddress('joe@example.net', 'Joe User');     // Add a recipient
$mail->addAddress('ellen@example.com');               // Name is optional
$mail->addReplyTo('info@example.com', 'Information');
$mail->addCC('cc@example.com');
$mail->addBCC('bcc@example.com');

$mail->addAttachment('/var/tmp/file.tar.gz');         // Add attachments
$mail->addAttachment('/tmp/image.jpg', 'new.jpg');    // Optional name
$mail->isHTML(true);                                  // Set email format to HTML

$mail->Subject = 'Here is the subject';
$mail->Body    = 'This is the HTML message body <b>in bold!</b>';
$mail->AltBody = 'This is the body in plain text for non-HTML mail clients';

if(!$mail->send()) {
    echo 'Message could not be sent.';
    echo 'Mailer Error: ' . $mail->ErrorInfo;
} else {
    echo 'Message has been sent';
}
```

You'll find plenty more to play with in the [examples](https://github.com/PHPMailer/PHPMailer/tree/master/examples) folder.

That's it. You should now be ready to use PHPMailer!

## Localization
PHPMailer defaults to English, but in the [language](https://github.com/PHPMailer/PHPMailer/tree/master/language/) folder you'll find numerous (46 at the time of writing!) translations for PHPMailer error messages that you may encounter. Their filenames contain [ISO 639-1](http://en.wikipedia.org/wiki/ISO_639-1) language code for the translations, for example `fr` for French. To specify a language, you need to tell PHPMailer which one to use, like this:

```php
// To load the French version
$mail->setLanguage('fr', '/optional/path/to/language/directory/');
```

We welcome corrections and new languages - if you're looking for corrections to do, run the [phpmailerLangTest.php](https://github.com/PHPMailer/PHPMailer/tree/master/test/phpmailerLangTest.php) script in the tests folder and it will show any missing translations.

## Documentation

Examples of how to use PHPMailer for common scenarios can be found in the [examples](https://github.com/PHPMailer/PHPMailer/tree/master/examples) folder. If you're looking for a good starting point, we recommend you start with [the Gmail example](https://github.com/PHPMailer/PHPMailer/tree/master/examples/gmail.phps).

There are tips and a troubleshooting guide in the [GitHub wiki](https://github.com/PHPMailer/PHPMailer/wiki). If you're having trouble, this should be the first place you look as it's the most frequently updated.

Complete generated API documentation is [available online](http://phpmailer.github.io/PHPMailer/).

You'll find some basic user-level docs in the [docs](docs/) folder, and you can generate complete API-level documentation using the [generatedocs.sh](https://github.com/PHPMailer/PHPMailer/tree/master/docs/generatedocs.sh) shell script in the docs folder, though you'll need to install [PHPDocumentor](http://www.phpdoc.org) first. You may find [the unit tests](https://github.com/PHPMailer/PHPMailer/tree/master/test/phpmailerTest.php) a good source of how to do various operations such as encryption.

If the documentation doesn't cover what you need, search the [many questions on Stack Overflow](http://stackoverflow.com/questions/tagged/phpmailer), and before you ask a question about "SMTP Error: Could not connect to SMTP host.", [read the troubleshooting guide](https://github.com/PHPMailer/PHPMailer/wiki/Troubleshooting).

## Tests

There is a PHPUnit test script in the [test](https://github.com/PHPMailer/PHPMailer/tree/master/test/) folder.

Build status: [![Build Status](https://travis-ci.org/PHPMailer/PHPMailer.svg)](https://travis-ci.org/PHPMailer/PHPMailer)

If this isn't passing, is there something you can do to help?

## Security

Please disclose any vulnerabilities found responsibly - report any security problems found to the maintainers privately.

PHPMailer versions prior to 5.2.14 (released November 2015) are vulnerable to [CVE-2015-8476](https://web.nvd.nist.gov/view/vuln/detail?vulnId=) an SMTP CRLF injection bug permitting arbitrary message sending.

PHPMailer versions prior to 5.2.10 (released May 2015) are vulnerable to [CVE-2008-5619](https://web.nvd.nist.gov/view/vuln/detail?vulnId=CVE-2008-5619), a remote code execution vulnerability in the bundled html2text library. This file was removed in 5.2.10, so if you are using a version prior to that and make use of the html2text function, it's vitally important that you upgrade and remove this file.

See [SECURITY](https://github.com/PHPMailer/PHPMailer/tree/master/SECURITY) for older security issues.

## Contributing

Please submit bug reports, suggestions and pull requests to the [GitHub issue tracker](https://github.com/PHPMailer/PHPMailer/issues).

We're particularly interested in fixing edge-cases, expanding test coverage and updating translations.

With the move to the PHPMailer GitHub organisation, you'll need to update any remote URLs referencing the old GitHub location with a command like this from within your clone:

```sh
git remote set-url upstream https://github.com/PHPMailer/PHPMailer.git
```

Please *don't* use the SourceForge or Google Code projects any more.

## Sponsorship

Development time and resources for PHPMailer are provided by [Smartmessages.net](https://info.smartmessages.net/), a powerful email marketing system.

<a href="https://info.smartmessages.net/"><img src="https://www.smartmessages.net/img/smartmessages-logo.svg" width="250" height="28" alt="Smartmessages email marketing"></a>

Other contributions are gladly received, whether in beer 🍺, T-shirts 👕, Amazon wishlist raids, or cold, hard cash 💰.

## Changelog

See [changelog](changelog.md).

## History
- PHPMailer was originally written in 2001 by Brent R. Matzelle as a [SourceForge project](http://sourceforge.net/projects/phpmailer/).
- Marcus Bointon (coolbru on SF) and Andy Prevost (codeworxtech) took over the project in 2004.
- Became an Apache incubator project on Google Code in 2010, managed by Jim Jagielski.
- Marcus created his fork on [GitHub](https://github.com/Synchro/PHPMailer).
- Jim and Marcus decide to join forces and use GitHub as the canonical and official repo for PHPMailer.
- PHPMailer moves to the [PHPMailer organisation](https://github.com/PHPMailer) on GitHub.

### What's changed since moving from SourceForge?
- Official successor to the SourceForge and Google Code projects.
- Test suite.
- Continuous integration with Travis-CI.
- Composer support.
- Public development.
- Additional languages and language strings.
- CRAM-MD5 authentication support.
- Preserves full repo history of authors, commits and branches from the original SourceForge project.
# ChangeLog

## Version 5.2.16 (June 6th 2016)
* Added DKIM example
* Fixed empty additional_parameters problem
* Fixed wrong version number in VERSION file!
* Improve line-length tests
* Use instance settings for SMTP::connect by default
* Use more secure auth mechanisms first

## Version 5.2.15 (May 10th 2016)
* Added ability to inject custom address validators, and set the default validator
* Fix TLS 1.2 compatibility
* Remove some excess line breaks in MIME structure
* Updated Polish, Russian, Brazilian Portuguese, Georgian translations
* More DRY!
* Improve error messages
* Update dependencies
* Add example showing how to handle multiple form file uploads
* Improve SMTP example
* Improve Windows compatibility
* Use consistent names for temp files
* Fix gmail XOAUTH2 scope, thanks to @sherryl4george
* Fix extra line break in getSentMIMEMessage()
* Improve DKIM signing to use SHA-2

## Version 5.2.14 (Nov 1st 2015)
* Allow addresses with IDN (Internationalized Domain Name) in PHP 5.3+, thanks to @fbonzon
* Allow access to POP3 errors
* Make all POP3 private properties and methods protected
* **SECURITY** Fix vulnerability that allowed email addresses with line breaks (valid in RFC5322) to pass to SMTP, permitting message injection at the SMTP level. Mitigated in both the address validator and in the lower-level SMTP class. Thanks to Takeshi Terada.
* Updated Brazilian Portuguese translations (Thanks to @phelipealves)

## Version 5.2.13 (Sep 14th 2015)
* Rename internal oauth class to avoid name clashes
* Improve Estonian translations

## Version 5.2.12 (Sep 1st 2015)
* Fix incorrect composer package dependencies
* Skip existing embedded image `cid`s in `msgHTML`

## Version 5.2.11 (Aug 31st 2015)
* Don't switch to quoted-printable for long lines if already using base64
* Fixed Travis-CI config when run on PHP 7
* Added Google XOAUTH2 authentication mechanism, thanks to @sherryl4george
* Add address parser for RFC822-format addresses
* Update MS Office MIME types
* Don't convert line breaks when using quoted-printable encoding
* Handle MS Exchange returning an invalid empty AUTH-type list in EHLO
* Don't set name or filename properties on MIME parts that don't have one

## Version 5.2.10 (May 4th 2015)
* Add custom header getter
* Use `application/javascript` for .js attachments
* Improve RFC2821 compliance for timelimits, especially for end-of-data
* Add Azerbaijani translations (Thanks to @mirjalal)
* Minor code cleanup for robustness
* Add Indonesian translations (Thanks to @ceceprawiro)
* Avoid `error_log` Debugoutput naming clash
* Add ability to parse server capabilities in response to EHLO (useful for SendGrid etc)
* Amended default values for WordWrap to match RFC
* Remove html2text converter class (has incompatible license)
* Provide new mechanism for injecting html to text converters
* Improve pointers to docs and support in README
* Add example file upload script
* Refactor and major cleanup of EasyPeasyICS, now a lot more usable
* Make set() method simpler and more reliable
* Add Malay translation (Thanks to @nawawi)
* Add Bulgarian translation (Thanks to @mialy)
* Add Armenian translation (Thanks to Hrayr Grigoryan)
* Add Slovenian translation (Thanks to Klemen Tušar)
* More efficient word wrapping
* Add support for S/MIME signing with additional CA certificate (thanks to @IgitBuh)
* Fix incorrect MIME structure when using S/MIME signing and isMail() (#372)
* Improved checks and error messages for missing extensions
* Store and report SMTP errors more consistently
* Add MIME multipart preamble for better Outlook compatibility
* Enable TLS encryption automatically if the server offers it
* Provide detailed errors when individual recipients fail
* Report more errors when connecting
* Add extras classes to composer classmap
* Expose stream_context_create options via new SMTPOptions property
* Automatic encoding switch to quoted-printable if message lines are too long
* Add Korean translation (Thanks to @ChalkPE)
* Provide a pointer to troubleshooting docs on SMTP connection failure

## Version 5.2.9 (Sept 25th 2014)
* **Important: The autoloader is no longer autoloaded by the PHPMailer class**
* Update html2text from https://github.com/mtibben/html2text
* Improve Arabic translations (Thanks to @tarekdj)
* Consistent handling of connection variables in SMTP and POP3
* PHPDoc cleanup
* Update composer to use PHPUnit 4.1
* Pass consistent params to callbacks
* More consistent handling of error states and debug output
* Use property defaults, remove constructors
* Remove unreachable code
* Use older regex validation pattern for troublesome PCRE library versions
* Improve PCRE detection in older PHP versions
* Handle debug output consistently, and always in UTF-8
* Allow user-defined debug output method via a callable
* msgHTML now converts data URIs to embedded images
* SMTP::getLastReply() will now always be populated
* Improved example code in README
* Ensure long filenames in Content-Disposition are encoded correctly
* Simplify SMTP debug output mechanism, clarify levels with constants
* Add SMTP connection check example
* Simplify examples, don't use mysql* functions

## Version 5.2.8 (May 14th 2014)
* Increase timeout to match RFC2821 section 4.5.3.2 and thus not fail greetdelays, fixes #104
* Add timestamps to default debug output
* Add connection events and new level 3 to debug output options
* Chinese language update (Thanks to @binaryoung)
* Allow custom Mailer types (Thanks to @michield)
* Cope with spaces around SMTP host specs
* Fix processing of multiple hosts in connect string
* Added Galician translation (Thanks to @donatorouco)
* Autoloader now prepends
* Docs updates
* Add Latvian translation (Thanks to @eddsstudio)
* Add Belarusian translation (Thanks to @amaksymiuk)
* Make autoloader work better on older PHP versions
* Avoid double-encoding if mbstring is overloading mail()
* Add Portuguese translation (Thanks to @Jonadabe)
* Make quoted-printable encoder respect line ending setting
* Improve Chinese translation (Thanks to @PeterDaveHello)
* Add Georgian translation (Thanks to @akalongman)
* Add Greek translation (Thanks to @lenasterg)
* Fix serverHostname on PHP < 5.3
* Improve performance of SMTP class
* Implement automatic 7bit downgrade
* Add Vietnamese translation (Thanks to @vinades)
* Improve example images, switch to PNG
* Add Croatian translation (Thanks to @hrvoj3e)
* Remove setting the Return-Path and deprecate the Return-path property - it's just wrong!
* Fix language file loading if CWD has changed (@stephandesouza)
* Add HTML5 email validation pattern
* Improve Turkish translations (Thanks to @yasinaydin)
* Improve Romanian translations (Thanks to @aflorea)
* Check php.ini for path to sendmail/qmail before using default
* Improve Farsi translation (Thanks to @MHM5000)
* Don't use quoted-printable encoding for multipart types
* Add Serbian translation (Thanks to ajevremovic at gmail.com)
* Remove useless PHP5 check
* Use SVG for build status badges
* Store MessageDate on creation
* Better default behaviour for validateAddress

## Version 5.2.7 (September 12th 2013)
* Add Ukrainian translation from @Krezalis
* Support for do_verp
* Fix bug in CRAM-MD5 AUTH
* Propagate Debugoutput option to SMTP class (@Reblutus)
* Determine MIME type of attachments automatically
* Add cross-platform, multibyte-safe pathinfo replacement (with tests) and use it
* Add a new 'html' Debugoutput type
* Clean up SMTP debug output, remove embedded HTML
* Some small changes in header formatting to improve IETF msglint test results
* Update test_script to use some recently changed features, rename to code_generator
* Generated code actually works!
* Update SyntaxHighlighter
* Major overhaul and cleanup of example code
* New PHPMailer graphic
* msgHTML now uses RFC2392-compliant content ids
* Add line break normalization function and use it in msgHTML
* Don't set unnecessary reply-to addresses
* Make fakesendmail.sh a bit cleaner and safer
* Set a content-transfer-encoding on multiparts (fixes msglint error)
* Fix cid generation in msgHTML (Thanks to @digitalthought)
* Fix handling of multiple SMTP servers (Thanks to @NanoCaiordo)
* SMTP->connect() now supports stream context options (Thanks to @stanislavdavid)
* Add support for iCal event alternatives (Thanks to @reblutus)
* Update to Polish language file (Thanks to Krzysztof Kowalewski)
* Update to Norwegian language file (Thanks to @datagutten)
* Update to Hungarian language file (Thanks to @dominicus-75)
* Add Persian/Farsi translation from @jaii
* Make SMTPDebug property type match type in SMTP class
* Add unit tests for DKIM
* Major refactor of SMTP class
* Reformat to PSR-2 coding standard
* Introduce autoloader
* Allow overriding of SMTP class
* Overhaul of PHPDocs
* Fix broken Q-encoding
* Czech language update (Thanks to @nemelu)
* Removal of excess blank lines in messages
* Added fake POP server and unit tests for POP-before-SMTP

## Version 5.2.6 (April 11th 2013)
* Reflect move to PHPMailer GitHub organisation at https://github.com/PHPMailer/PHPMailer
* Fix unbumped version numbers
* Update packagist.org with new location
* Clean up Changelog

## Version 5.2.5 (April 6th 2013)
* First official release after move from Google Code
* Fixes for qmail when sending via mail()
* Merge in changes from Google code 5.2.4 release
* Minor coding standards cleanup in SMTP class
* Improved unit tests, now tests S/MIME signing
* Travis-CI support on GitHub, runs tests with fake SMTP server

## Version 5.2.4 (February 19, 2013)
* Fix tag and version bug.
* un-deprecate isSMTP(), isMail(), IsSendmail() and isQmail().
* Numerous translation updates

## Version 5.2.3 (February 8, 2013)
* Fix issue with older PCREs and ValidateAddress() (Bugz: 124)
* Add CRAM-MD5 authentication, thanks to Elijah madden, https://github.com/okonomiyaki3000
* Replacement of obsolete Quoted-Printable encoder with a much better implementation
* Composer package definition
* New language added: Hebrew

## Version 5.2.2 (December 3, 2012)
* Some fixes and syncs from https://github.com/Synchro/PHPMailer
* Add Slovak translation, thanks to Michal Tinka

## Version 5.2.2-rc2 (November 6, 2012)
* Fix SMTP server rotation (Bugz: 118)
* Allow override of autogen'ed 'Date' header (for Drupal's
  og_mailinglist module)
* No whitespace after '-f' option (Bugz: 116)
* Work around potential warning (Bugz: 114)

## Version 5.2.2-rc1 (September 28, 2012)
* Header encoding works with long lines (Bugz: 93)
* Turkish language update (Bugz: 94)
* undefined $pattern in EncodeQ bug squashed (Bugz: 98)
* use of mail() in safe_mode now works (Bugz: 96)
* ValidateAddress() now 'public static' so people can override the
  default and use their own validation scheme.
* ValidateAddress() no longer uses broken FILTER_VALIDATE_EMAIL
* Added in AUTH PLAIN SMTP authentication

## Version 5.2.2-beta2 (August 17, 2012)
* Fixed Postfix VERP support (Bugz: 92)
* Allow action_function callbacks to pass/use
  the From address (passed as final param)
* Prevent inf look for get_lines() (Bugz: 77)
* New public var ($UseSendmailOptions). Only pass sendmail()
  options iff we really are using sendmail or something sendmail
  compatible. (Bugz: 75)
* default setting for LE returned to "\n" due to popular demand.

## Version 5.2.2-beta1 (July 13, 2012)
* Expose PreSend() and PostSend() as public methods to allow
  for more control if serializing message sending.
* GetSentMIMEMessage() only constructs the message copy when
 needed. Save memory.
* Only pass params to mail() if the underlying MTA is
  "sendmail" (as defined as "having the string sendmail
  in its pathname") [#69]
* Attachments now work with Amazon SES and others [Bugz#70]
* Debug output now sent to stdout (via echo) or error_log [Bugz#5]
* New var: Debugoutput (for above) [Bugz#5]
* SMTP reads now Timeout aware (new var: Timeout=15) [Bugz#71]
* SMTP reads now can have a Timelimit associated with them
  (new var: Timelimit=30)[Bugz#71]
* Fix quoting issue associated with charsets
* default setting for LE is now RFC compliant: "\r\n"
* Return-Path can now be user defined (new var: ReturnPath)
  (the default is "" which implies no change from previous
  behavior, which was to use either From or Sender) [Bugz#46]
* X-Mailer header can now be disabled (by setting to a
  whitespace string, eg "  ") [Bugz#66]
* Bugz closed: #68, #60, #42, #43, #59, #55, #66, #48, #49,
               #52, #31, #41, #5. #70, #69

## Version 5.2.1 (January 16, 2012)
* Closed several bugs #5
* Performance improvements
* MsgHTML() now returns the message as required.
* New method: GetSentMIMEMessage() (returns full copy of sent message)

## Version 5.2 (July 19, 2011)
* protected MIME body and header
* better DKIM DNS Resource Record support
* better aly handling
* htmlfilter class added to extras
* moved to Apache Extras

## Version 5.1 (October 20, 2009)
* fixed filename issue with AddStringAttachment (thanks to Tony)
* fixed "SingleTo" property, now works with Senmail, Qmail, and SMTP in
  addition to PHP mail()
* added DKIM digital signing functionality, new properties:
  - DKIM_domain (sets the domain name)
  - DKIM_private (holds DKIM private key)
  - DKIM_passphrase (holds your DKIM passphrase)
  - DKIM_selector (holds the DKIM "selector")
  - DKIM_identity (holds the identifying email address)
* added callback function support
  - callback function parameters include:
    result, to, cc, bcc, subject and body
  - see the test/test_callback.php file for usage.
* added "auto" identity functionality
  - can automatically add:
    - Return-path (if Sender not set)
    - Reply-To (if ReplyTo not set)
  - can be disabled:
    - $mail->SetFrom('yourname@yourdomain.com','First Last',false);
    - or by adding the $mail->Sender and/or $mail->ReplyTo properties

Note: "auto" identity added to help with emails ending up in spam or junk boxes because of missing headers

## Version 5.0.2 (May 24, 2009)
* Fix for missing attachments when inline graphics are present
* Fix for missing Cc in header when using SMTP (mail was sent,
  but not displayed in header -- Cc receiver only saw email To:
  line and no Cc line, but did get the email (To receiver
  saw same)

## Version 5.0.1 (April 05, 2009)
* Temporary fix for missing attachments

## Version 5.0.0 (April 02, 2009)
With the release of this version, we are initiating a new version numbering
system to differentiate from the PHP4 version of PHPMailer.
Most notable in this release is fully object oriented code.

### class.smtp.php:
* Refactored class.smtp.php to support new exception handling
* code size reduced from 29.2 Kb to 25.6 Kb
* Removed unnecessary functions from class.smtp.php:
  - public function Expand($name) {
  - public function Help($keyword="") {
  - public function Noop() {
  - public function Send($from) {
  - public function SendOrMail($from) {
  - public function Verify($name) {

###  class.phpmailer.php:
* Refactored class.phpmailer.php with new exception handling
* Changed processing functionality of Sendmail and Qmail so they cannot be
  inadvertently used
* removed getFile() function, just became a simple wrapper for
  file_get_contents()
* added check for PHP version (will gracefully exit if not at least PHP 5.0)
* enhanced code to check if an attachment source is the same as an embedded or
  inline graphic source to eliminate duplicate attachments

### New /test_script
We have written a test script you can use to test the script as part of your
installation. Once you press submit, the test script will send a multi-mime
email with either the message you type in or an HTML email with an inline
graphic. Two attachments are included in the email (one of the attachments
is also the inline graphic so you can see that only one copy of the graphic
is sent in the email). The test script will also display the functional
script that you can copy/paste to your editor to duplicate the functionality.

### New examples
All new examples in both basic and advanced modes. Advanced examples show
   Exception handling.

### PHPDocumentator (phpdocs) documentation for PHPMailer version 5.0.0
All new documentation

## Version 2.3 (November 06, 2008)
* added Arabic language (many thanks to Bahjat Al Mostafa)
* removed English language from language files and made it a default within
  class.phpmailer.php - if no language is found, it will default to use
  the english language translation
* fixed public/private declarations
* corrected line 1728, $basedir to $directory
* added $sign_cert_file to avoid improper duplicate use of $sign_key_file
* corrected $this->Hello on line 612 to $this->Helo
* changed default of $LE to "\r\n" to comply with RFC 2822. Can be set by the user
  if default is not acceptable
* removed trim() from return results in EncodeQP
* /test and three files it contained are removed from version 2.3
* fixed phpunit.php for compliance with PHP5
* changed $this->AltBody = $textMsg; to $this->AltBody = html_entity_decode($textMsg);
* We have removed the /phpdoc from the downloads. All documentation is now on
  the http://phpmailer.codeworxtech.com website.

## Version 2.2.1 () July 19 2008
* fixed line 1092 in class.smtp.php (my apologies, error on my part)

## Version 2.2 () July 15 2008
* Fixed redirect issue (display of UTF-8 in thank you redirect)
* fixed error in getResponse function declaration (class.pop3.php)
* PHPMailer now PHP6 compliant
* fixed line 1092 in class.smtp.php (endless loop from missing = sign)

## Version 2.1 (Wed, June 04 2008)
NOTE: WE HAVE A NEW LANGUAGE VARIABLE FOR DIGITALLY SIGNED S/MIME EMAILS. IF YOU CAN HELP WITH LANGUAGES OTHER THAN ENGLISH AND SPANISH, IT WOULD BE APPRECIATED.

* added S/MIME functionality (ability to digitally sign emails)
  BIG THANKS TO "sergiocambra" for posting this patch back in November 2007.
  The "Signed Emails" functionality adds the Sign method to pass the private key
  filename and the password to read it, and then email will be sent with
  content-type multipart/signed and with the digital signature attached.
* fully compatible with E_STRICT error level
  - Please note:
    In about half the test environments this development version was subjected
    to, an error was thrown for the date() functions used (line 1565 and 1569).
    This is NOT a PHPMailer error, it is the result of an incorrectly configured
    PHP5 installation. The fix is to modify your 'php.ini' file and include the
    date.timezone = Etc/UTC (or your own zone)
    directive, to your own server timezone
  - If you do get this error, and are unable to access your php.ini file:
    In your PHP script, add
    `date_default_timezone_set('Etc/UTC');`
  - do not try to use
    `$myVar = date_default_timezone_get();`
    as a test, it will throw an error.
* added ability to define path (mainly for embedded images)
  function `MsgHTML($message,$basedir='')` ... where:
  `$basedir` is the fully qualified path
* fixed `MsgHTML()` function:
  - Embedded Images where images are specified by `<protocol>://` will not be altered or embedded
* fixed the return value of SMTP exit code ( pclose )
* addressed issue of multibyte characters in subject line and truncating
* added ability to have user specified Message ID
  (default is still that PHPMailer create a unique Message ID)
* corrected unidentified message type to 'application/octet-stream'
* fixed chunk_split() multibyte issue (thanks to Colin Brown, et al).
* added check for added attachments
* enhanced conversion of HTML to text in MsgHTML (thanks to "brunny")

## Version 2.1.0beta2 (Sun, Dec 02 2007)
* implemented updated EncodeQP (thanks to coolbru, aka Marcus Bointon)
* finished all testing, all known bugs corrected, enhancements tested

Note: will NOT work with PHP4.

Please note, this is BETA software **DO NOT USE THIS IN PRODUCTION OR LIVE PROJECTS; INTENDED STRICTLY FOR TESTING**

## Version 2.1.0beta1
Please note, this is BETA software
** DO NOT USE THIS IN PRODUCTION OR LIVE PROJECTS
 INTENDED STRICTLY FOR TESTING

## Version 2.0.0 rc2 (Fri, Nov 16 2007), interim release
* implements new property to control VERP in class.smtp.php
  example (requires instantiating class.smtp.php):
  $mail->do_verp = true;
* POP-before-SMTP functionality included, thanks to Richard Davey
  (see class.pop3.php & pop3_before_smtp_test.php for examples)
* included example showing how to use PHPMailer with GMAIL
* fixed the missing Cc in SendMail() and Mail()

## Version 2.0.0 rc1 (Thu, Nov 08 2007), interim release
* dramatically simplified using inline graphics ... it's fully automated and requires no user input
* added automatic document type detection for attachments and pictures
* added MsgHTML() function to replace Body tag for HTML emails
* fixed the SendMail security issues (input validation vulnerability)
* enhanced the AddAddresses functionality so that the "Name" portion is used in the email address
* removed the need to use the AltBody method (set from the HTML, or default text used)
* set the PHP Mail() function as the default (still support SendMail, SMTP Mail)
* removed the need to set the IsHTML property (set automatically)
* added Estonian language file by Indrek P&auml;ri
* added header injection patch
* added "set" method to permit users to create their own pseudo-properties like 'X-Headers', etc.
* fixed warning message in SMTP get_lines method
* added TLS/SSL SMTP support.
* PHPMailer has been tested with PHP4 (4.4.7) and PHP5 (5.2.7)
* Works with PHP installed as a module or as CGI-PHP
NOTE: will NOT work with PHP5 in E_STRICT error mode

## Version 1.73 (Sun, Jun 10 2005)
* Fixed denial of service bug: http://www.cybsec.com/vuln/PHPMailer-DOS.pdf
* Now has a total of 20 translations
* Fixed alt attachments bug: http://tinyurl.com/98u9k

## Version 1.72 (Wed, May 25 2004)
* Added Dutch, Swedish, Czech, Norwegian, and Turkish translations.
* Received: Removed this method because spam filter programs like
  SpamAssassin reject this header.
* Fixed error count bug.
* SetLanguage default is now "language/".
* Fixed magic_quotes_runtime bug.

## Version 1.71 (Tue, Jul 28 2003)
* Made several speed enhancements
* Added German and Italian translation files
* Fixed HELO/AUTH bugs on keep-alive connects
* Now provides an error message if language file does not load
* Fixed attachment EOL bug
* Updated some unclear documentation
* Added additional tests and improved others

## Version 1.70 (Mon, Jun 20 2003)
* Added SMTP keep-alive support
* Added IsError method for error detection
* Added error message translation support (SetLanguage)
* Refactored many methods to increase library performance
* Hello now sends the newer EHLO message before HELO as per RFC 2821
* Removed the boundary class and replaced it with GetBoundary
* Removed queue support methods
* New $Hostname variable
* New Message-ID header
* Received header reformat
* Helo variable default changed to $Hostname
* Removed extra spaces in Content-Type definition (#667182)
* Return-Path should be set to Sender when set
* Adds Q or B encoding to headers when necessary
* quoted-encoding should now encode NULs \000
* Fixed encoding of body/AltBody (#553370)
* Adds "To: undisclosed-recipients:;" when all recipients are hidden (BCC)
* Multiple bug fixes

## Version 1.65 (Fri, Aug 09 2002)
* Fixed non-visible attachment bug (#585097) for Outlook
* SMTP connections are now closed after each transaction
* Fixed SMTP::Expand return value
* Converted SMTP class documentation to phpDocumentor format

## Version 1.62 (Wed, Jun 26 2002)
* Fixed multi-attach bug
* Set proper word wrapping
* Reduced memory use with attachments
* Added more debugging
* Changed documentation to phpDocumentor format

## Version 1.60 (Sat, Mar 30 2002)
* Sendmail pipe and address patch (Christian Holtje)
* Added embedded image and read confirmation support (A. Ognio)
* Added unit tests
* Added SMTP timeout support (*nix only)
* Added possibly temporary PluginDir variable for SMTP class
* Added LE message line ending variable
* Refactored boundary and attachment code
* Eliminated SMTP class warnings
* Added SendToQueue method for future queuing support

## Version 1.54 (Wed, Dec 19 2001)
* Add some queuing support code
* Fixed a pesky multi/alt bug
* Messages are no longer forced to have "To" addresses

## Version 1.50 (Thu, Nov 08 2001)
* Fix extra lines when not using SMTP mailer
* Set WordWrap variable to int with a zero default

## Version 1.47 (Tue, Oct 16 2001)
* Fixed Received header code format
* Fixed AltBody order error
* Fixed alternate port warning

## Version 1.45 (Tue, Sep 25 2001)
* Added enhanced SMTP debug support
* Added support for multiple ports on SMTP
* Added Received header for tracing
* Fixed AddStringAttachment encoding
* Fixed possible header name quote bug
* Fixed wordwrap() trim bug
* Couple other small bug fixes

## Version 1.41 (Wed, Aug 22 2001)
* Fixed AltBody bug w/o attachments
* Fixed rfc_date() for certain mail servers

## Version 1.40 (Sun, Aug 12 2001)
* Added multipart/alternative support (AltBody)
* Documentation update
* Fixed bug in Mercury MTA

## Version 1.29 (Fri, Aug 03 2001)
* Added AddStringAttachment() method
* Added SMTP authentication support

## Version 1.28 (Mon, Jul 30 2001)
* Fixed a typo in SMTP class
* Fixed header issue with Imail (win32) SMTP server
* Made fopen() calls for attachments use "rb" to fix win32 error

## Version 1.25 (Mon, Jul 02 2001)
* Added RFC 822 date fix (Patrice)
* Added improved error handling by adding a $ErrorInfo variable
* Removed MailerDebug variable (obsolete with new error handler)

## Version 1.20 (Mon, Jun 25 2001)
* Added quoted-printable encoding (Patrice)
* Set Version as public and removed PrintVersion()
* Changed phpdoc to only display public variables and methods

## Version 1.19 (Thu, Jun 21 2001)
* Fixed MS Mail header bug
* Added fix for Bcc problem with mail(). *Does not work on Win32*
  (See PHP bug report: http://www.php.net/bugs.php?id=11616)
* mail() no longer passes a fifth parameter when not needed

## Version 1.15 (Fri, Jun 15 2001)
Note: these changes contributed by Patrice Fournier
* Changed all remaining \n to \r\n
* Bcc: header no longer written to message except
  when sent directly to sendmail
* Added a small message to non-MIME compliant mail reader
* Added Sender variable to change the Sender email
  used in -f for sendmail/mail and in 'MAIL FROM' for smtp mode
* Changed boundary setting to a place it will be set only once
* Removed transfer encoding for whole message when using multipart
* Message body now uses Encoding in multipart messages
* Can set encoding and type to attachments 7bit, 8bit
  and binary attachment are sent as is, base64 are encoded
* Can set Encoding to base64 to send 8 bits body
  through 7 bits servers

## Version 1.10 (Tue, Jun 12 2001)
* Fixed win32 mail header bug (printed out headers in message body)

## Version 1.09 (Fri, Jun 08 2001)
* Changed date header to work with Netscape mail programs
* Altered phpdoc documentation

## Version 1.08 (Tue, Jun 05 2001)
* Added enhanced error-checking
* Added phpdoc documentation to source

## Version 1.06 (Fri, Jun 01 2001)
* Added optional name for file attachments

## Version 1.05 (Tue, May 29 2001)
* Code cleanup
* Eliminated sendmail header warning message
* Fixed possible SMTP error

## Version 1.03 (Thu, May 24 2001)
* Fixed problem where qmail sends out duplicate messages

## Version 1.02 (Wed, May 23 2001)
* Added multiple recipient and attachment Clear* methods
* Added Sendmail public variable
* Fixed problem with loading SMTP library multiple times

## Version 0.98 (Tue, May 22 2001)
* Fixed problem with redundant mail hosts sending out multiple messages
* Added additional error handler code
* Added AddCustomHeader() function
* Added support for Microsoft mail client headers (affects priority)
* Fixed small bug with Mailer variable
* Added PrintVersion() function

## Version 0.92 (Tue, May 15 2001)
* Changed file names to class.phpmailer.php and class.smtp.php to match
  current PHP class trend.
* Fixed problem where body not being printed when a message is attached
* Several small bug fixes

## Version 0.90 (Tue, April 17 2001)
* Initial public release
#PHPMailer Extras

These classes provide optional additional functions to PHPMailer.

These are not loaded by the PHPMailer autoloader, so in some cases you may need to `require` them yourself before using them.

##EasyPeasyICS

This class was originally written by Manuel Reinhard and provides a simple means of generating ICS/vCal files that are used in sending calendar events. PHPMailer does not use it directly, but you can use it to generate content appropriate for placing in the `Ical` property of PHPMailer. The PHPMailer project is now its official home as Manuel has given permission for that and is no longer maintaining it himself.

##htmlfilter

This class by Konstantin Riabitsev and Jim Jagielski implements HTML filtering to remove potentially malicious tags, such as `<script>` or `onclick=` attributes that can result in XSS attacks. This is a simple filter and is not as comprehensive as [HTMLawed](http://www.bioinformatics.org/phplabware/internal_utilities/htmLawed/) or [HTMLPurifier](http://htmlpurifier.org), but it's easier to use and considerably better than nothing! PHPMailer does not use it directly, but you may want to apply it to user-supplied HTML before using it as a message body.

##NTLM_SASL_client

This class by Manuel Lemos (bundled with permission) adds the ability to authenticate with Microsoft Windows mail servers that use NTLM-based authentication. It is used by PHPMailer if you send via SMTP and set the `AuthType` property to `NTLM`; you will also need to use the `Realm` and `Workstation` properties. The original source is [here](http://www.phpclasses.org/browse/file/7495.html).
datetimepicker
==============
[Documentation][doc]


jQuery Plugin Date and Time Picker

DateTimePicker

![ScreenShot](https://raw.github.com/xdan/datetimepicker/master/screen/1.png)

DatePicker

![ScreenShot](https://raw.github.com/xdan/datetimepicker/master/screen/2.png)

TimePicker

![ScreenShot](https://raw.github.com/xdan/datetimepicker/master/screen/3.png)

Options to highlight individual dates or periods

![ScreenShot](https://raw.github.com/Mingpao/datetimepicker/master/screen/4.png)

![ScreenShot](https://raw.github.com/Mingpao/datetimepicker/master/screen/5.png)

![ScreenShot](https://raw.github.com/Mingpao/datetimepicker/master/screen/6.png)

[doc]: http://xdsoft.net/jqplugins/datetimepicker/
This is where language files should be placed.

Please DO NOT translate these directly use this service: https://www.transifex.com/projects/p/tinymce/
lasso
=========

lasso.js is a D3 plugin that allows you to tag elements on a page by drawing a line over or around objects. Functions can be run based on the lasso action. This functionality can be useful for brushing or filtering.

An example of the lasso implemented in a scatterplot can be found here: [http://bl.ocks.org/skokenes/511c5b658c405ad68941](http://bl.ocks.org/skokenes/511c5b658c405ad68941)

This example is based off of Mike Bostock's scatterplot example here: [http://bl.ocks.org/mbostock/3887118](http://bl.ocks.org/mbostock/3887118)

Lassoing tags
--
When the lasso is used, it tags elements by adding properties to their data. The properties are:

- possible: while drawing a lasso, if an element is part of the final selection that would be made if the lasso was completed at that instance, this value is true. Otherwise, it is false.
- selected: when a lasso is completed, all elements that were tagged as possible are given a selected value of true. Otherwise, the value is false.

The tags can be used in combination with functions to perform actions like styling the possible or selected values while the lasso is in use.

Note that the lasso only works with elements whose data is defined as an object.


Function Overview
--
**d3.lasso**()

Creates a new lasso object. This object can then have parameters set before the lasso is drawn.
```
var lasso = d3.lasso(); // creates a new lasso
```

lasso.**items**(_[selection]_)

The items() parameter takes in a d3 selection. Each element in the selection will be tagged with lasso-specific properties when the lasso is used. If no input is specified, the function returns the lasso's current items.
```
lasso.items(d3.selectAll("circle")); // sets all circles on the page to be lasso-able
```

lasso.**hoverSelect**(_[bool]_)

The hoverSelect() parameter takes in a boolean that determines whether objects can be lassoed by hovering over an element during lassoing. The default value is set to true. If no input is specified, the function returns the lasso's current hover parameter.
```
lasso.hoverSelect(true); // allows hovering of elements for selection during lassoing
```

lasso.**closePathSelect**(_[bool]_)

The closePathSelect() parameter takes in a boolean that determines whether objects can be lassoed by drawing a loop around them. The default value is set to true. If no input is specified, the function returns the lasso's current parameter.
```
lasso.closePathSelect(true); // allows looping of elements for selection during lassoing
```

lasso.**closePathDistance**(_[num]_)

The closePathDistance() parameter takes in a number that specifies the maximum distance in pixels from the lasso origin that a lasso needs to be drawn in order to complete the loop and select elements. This parameter only works if closePathSelect is set to true; If no input is specified, the function returns the lasso's current parameter.
```
lasso.closePathDistance(75); // the lasso loop will complete itself whenever the lasso end is within 75 pixels of the origin
```

lasso.**area**(_[sel]_)

The area() parameter takes in a selection representing the element to be used as a target area for the lasso event. If no input is specified, the function returns the current area selection.
```
lasso.area(d3.select("#myLassoRect")); // the lasso will be trigger whenever a user clicks and drags on #myLassoRect
```

lasso.**on**(_type,[func]_)

The on() parameter takes in a type of event and a function for that event. There are 3 types of events that can be defined:
- start: this function will be executed whenever a lasso is started
- draw: this function will execute repeatedly as the lasso is drawn
- end: this function will be executed whenever a lasso is completed

If no function is specified, the function will return the current function defined for the type specified.
```
lasso.on("start",function() { alert("lasso started!"); }); // every time a lasso is started, an alert will trigger
```

Initiating a lasso
--
Once a lasso object is defined, it can be added to a page by calling it on an element like an svg.
```
var lasso = d3.lasso()
                .items(d3.selectAll("circle")) // Create a lasso and provide it some target elements
                .area(de.select("#myLassoRect")); // Sets the drag area for the lasso on the rectangle #myLassoRect
d3.select("svg").call(lasso); // Initiate the lasso on an svg element
```

If a lasso is going to be used on graphical elements that have been translated via a g element acting as a container, which is a common practice for incorporating chart margins, then the lasso should be called on that g element so that it is in the same coordinate system as the graphical elements.# DataTables plug-in for jQuery

DataTables is a table enhancing plug-in for the [jQuery](//jquery.com) Javascript library, adding sorting, paging and filtering abilities to plain HTML tables with minimal effort. The stated goal of DataTables is:

> To enhance the accessibility of data in HTML tables.

To meet this goal, DataTables is developed with two distinct groups of users in mind:

* You the developers using DataTables. For developers DataTables provides a wide array of options for how data should be obtained, displayed and acted upon, along with an extensive API for accessing and manipulating the table.

* End users. For those using the interface DataTables presents, actions to get the most from the information contained in tables, such as sorting and filtering, along with paging and scrolling of the data in table, are easy to use, intuitive and fast.


## Installing DataTables

To use DataTables, the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](//datatables.net/manual/installation) for full details.

### NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


## Usage

In its simplest case, DataTables can be initialised with a single line of Javascript:

```js
$('table').dataTable();
```

where the jQuery selector is used to obtain a reference to the table you want to enhance with DataTables. Optional configuration parameters can be passed in to DataTables to have it perform certain actions by using a configuration object as the parameter passed in to the DataTables constructor. For example:

```js
$('table').dataTable( {
  paginate: false,
  scrollY: 300
} );
```

will disable paging and enable scrolling.

A full list of the options available for DataTables are available in the [documentation](//datatables.net).


## Documentation

Full documentation of the DataTables options, API and plug-in interface are available on the [DataTables web-site](//datatables.net). The site also contains information on the wide variety of plug-ins that are available for DataTables, which can be used to enhance and customise your table even further.


## Support

Support for DataTables is available through the [DataTables forums](//datatables.net/forums) and [commercial support options](//datatables.net/support) are available.


## License

DataTables is release under the [MIT license](//datatables.net/license). You are free to use, modify and distribute this software, as long as the copyright header is left intact (specifically the comment block which starts with `/*!`.
Please post support requests and questions in the DataTables forums at https://datatables.net/forums. Support requests posted here will be closed. This allows all questions to be located in a single, searchable, location.

When you post your question in the DataTables forums, please ensure that you include a link to the page showing the issue so it can be debugged.
# Support requests

Please direct support requests to the [DataTables forums](https://datatables.net/forums), ensuring that you provide a link to a test page that shows the problem and a full description of the issue. If you require urgent help, [priority support](https://datatables.net/support) is available.


# Contributing code

If you are thinking of contributing code to DataTables, first of all, thank you! All fixes, patches and enhancements to DataTables are very warmly welcomed. In order to keep thing manageable, there are a number of guidelines that should be followed in order to ensure that your modification is included in DataTables as quickly as possible:

1. Make contributions in the DataTables/DataTablesSrc repo. Changes to the built files in the built repo (DataTables/DataTables) will not be accepted since they would be overwritten by the next build!

2. Follow the style of the code in the existing files. They might not be to everyone's tastes, but consistency is key for a mature project like DataTables. DataTables doesn't have a coding standards document, but simple common sense of following the same style as in the existing files is ideal. For example use tabs not spaces (as you will see all source files use tabs).

3. Link to a test page showing the bug you are fixing or the feature you are adding. This allows to me to quickly identify what is being changed and why. Don't worry about being verbose in pull requests - its much better to know exactly what is changing and why!

4. DataTables is a large and complex project and it isn't always possible or suitable to pull in every suggested change. Please don't be offended if a pull request is not merged in, it will explained why not if this is the case. Also it isn't always possible to fully check and test pull requests as quickly as I would like due to other commitments. Again this is no reflection on your pull request, just the busy life that we all lead! If you have any questions about your potential contribution and its place in the DataTables project structure, please ask ahead of time in the [DataTables forums](//datatables.net/forums).

5. Pull requests will only be accepted if you acknowledge that your contribution is offered under and will be made available under the project's existing license (MIT). If your initial pull request doesn't explicitly acknowledge this I'll ask before it is pulled in.# KeyTable

KeyTable provides Excel like cell navigation on any table. Events (focus, blur, action etc) can be assigned to individual cells, columns, rows or all cells.


# Installation

To use KeyTable the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/keytable/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-keytable`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

KeyTable is initialised using the `keys` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details.

```js
$(document).ready( function () {
    $('#myTable').DataTable( {
    	keys: true
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/keytable/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of KeyTable and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/KeyTable).

# FixedColumns

When making use of DataTables' x-axis scrolling feature (`scrollX`), you may wish to fix the left or right most columns in place. This extension for DataTables provides exactly this option (for non-scrolling tables, please use the FixedHeader extension, which can fix the header and footer).


# Installation

To use FixedColumns the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/fixedcolumns/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-fixedcolumns`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

FixedColumns is initialised using the `fixedColumns` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details. DataTables' scrolling options should also be enabled to use this feature.

Example:

```js
$(document).ready(function() {
	var table = $('#example').DataTable( {
		scrollY:        "300px",
		scrollX:        true,
		scrollCollapse: true,
		paging:         false,
		fixedColumns:   true
	} );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/fixedcolumns/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of FixedColumns and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/FixedColumns).
# FixedHeader

The FixedHeader plug-in will freeze in place the header, footer and left and/or right most columns in a DataTable, ensuring that title information will remain always visible.


# Installation

To use FixedHeader the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/fixedheader/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-fixedheader`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

FixedHeader is initialised using the `fixedHeader` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details.

Example:

```js
$(document).ready( function () {
    $('#myTable').DataTable( {
    	fixedHeader: true
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/fixedheader/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of FixedHeader and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/FixedHeader).

# Buttons

The Buttons extension for DataTables provides a common set of options, API methods and styling to display buttons on a page that will interact with a DataTable. Modules are also provided for data export, printing and column visibility control.


# Installation

To use Buttons the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/buttons/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-buttons`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

Buttons is initialised using the `buttons` option in the DataTables constructor, giving an array of the buttons that should be shown. Further options can be specified using this option as an object - see the documentation for details. For example:

```js
$(document).ready( function () {
    $('#example').DataTable( {
    	buttons: [ 'csv', 'excel', 'pdf', 'print' ]
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/buttons/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of Buttons and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/Buttons)

# AutoFill

AutoFill adds an Excel data fill like option to a DataTable to click and drag over multiple cells, filling in information over the selected cells and incrementing numbers as needed.


# Installation

To use AutoFill the best way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/autofill/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-autofill`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

AutoFill is initialised using the `autoFill` option in the DataTables constructor. Further options can be specified using this option as an object - see the documentation for details. For example:

```js
$(document).ready( function () {
    $('#example').DataTable( {
    	autoFill: true
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/autofill/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of AutoFill and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/AutoFill)

# Select

Select adds item selection capabilities to a DataTable. Items can be rows, columns or cells, which can be selected independently, or together. Item selection can be particularly useful in interactive tables where users can perform some action on the table such as editing.


# Installation

To use Select the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/select/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-select`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

Select is initialised using the `select` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details.

Example:

```js
$(document).ready( function () {
    $('#myTable').DataTable( {
    	select: true
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/select/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of Select and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/Select).

# RowReorder

RowReorder adds the ability for rows in a DataTable to be reordered through user interaction with the table (click and drag / touch and drag). Integration with Editor's multi-row editing feature is also available to update rows immediately. 


# Installation

To use RowReorder the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/rowreorder/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-rowreorder`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

RowReorder is initialised using the `rowReorder` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details.

Example:

```js
$(document).ready( function () {
    $('#myTable').DataTable( {
    	rowReorder: true
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/rowreorder/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of RowReorder and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/RowReorder).

# ColReorder

ColReorder adds the ability for the end user to click and drag column headers to reorder a table as they see fit, to DataTables. See the [documentation](http://datatables.net/extensions/colreorder/) for full details.


# Installation

To use ColReorder the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/colreorder/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-colreorder`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

ColReorder is initialised using the `colReorder` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details.

Example:

```js
$(document).ready( function () {
    $('#myTable').DataTable( {
    	colReorder: true
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/colreorder/)
* [DataTables support forums](http://datatables.net/forums)
# Scroller

Scroller is a virtual rendering plug-in for DataTables which allows large datasets to be drawn on screen every quickly. What the virtual rendering means is that only the visible portion of the table (and a bit to either side to make the scrolling smooth) is drawn, while the scrolling container gives the visual impression that the whole table is visible. This is done by making use of the pagination abilities of DataTables and moving the table around in the scrolling container DataTables adds to the page. The scrolling container is forced to the height it would be for the full table display using an extra element.

Key features include:

* Speed! The aim of Scroller for DataTables is to make rendering large data sets fast
* Full compatibility with DataTables' deferred rendering for maximum speed
* Integration with state saving in DataTables (scrolling position is saved)
* Support for scrolling with millions of rows
* Easy to use


# Installation

To use Scroller the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/scroller/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-scroller`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

Scroller is initialised using the `scroller` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details.

```js
$(document).ready( function () {
	$('#example').DataTable( {
		scroller: true
	} );
} );
```

Note that rows in the table must all be the same height. Information in a cell which expands on to multiple lines will cause some odd behaviour in the scrolling. Additionally, the table's `cellspacing` parameter must be set to 0, again to ensure the information display is correct.


# Documentation / support

* [Documentation](https://datatables.net/extensions/scroller/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of Scroller and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/Scroller)

# Responsive

Responsive will automatically optimise the table's layout for different screen sizes through the dynamic column visibility control, making your tables useful on desktop and mobile screens.


# Installation

To use Responsive the primary way to obtain the software is to use the [DataTables downloader](//datatables.net/download). You can also include the individual files from the [DataTables CDN](//cdn.datatables.net). See the [documentation](http://datatables.net/extensions/responsive/) for full details.

## NPM and Bower

If you prefer to use a package manager such as NPM or Bower, distribution repositories are available with software built from this repository under the name `datatables.net-responsive`. Styling packages for Bootstrap, Foundation and other styling libraries are also available by adding a suffix to the package name.

Please see the DataTables [NPM](//datatables.net/download/npm) and [Bower](//datatables.net/download/bower) installation pages for further information. The [DataTables installation manual](//datatables.net/manual/installation) also has details on how to use package managers with DataTables.


# Basic usage

Responsive is initialised using the `responsive` option in the DataTables constructor - a simple boolean `true` will enable the feature. Further options can be specified using this option as an object - see the documentation for details.

Example:

```js
$(document).ready( function () {
    $('#myTable').DataTable( {
    	responsive: true
    } );
} );
```


# Documentation / support

* [Documentation](https://datatables.net/extensions/responsive/)
* [DataTables support forums](http://datatables.net/forums)


# GitHub

If you fancy getting involved with the development of Responsive and help make it better, please refer to its [GitHub repo](https://github.com/DataTables/Responsive).

