# signature-commons-ui
A front-end UI for demoing API integration. Currently available at: http://amp.pharm.mssm.edu/sigcom/

## Development
Before starting, install the project dependencies:

```bash
npm install
```

### `npm run dev`

Runs the app in the development mode.
Open [http://localhost:3000](http://localhost:3000) to view it in the browser.

The page will reload if you make edits.<br>
You will also see any lint errors in the console.

### `npm run export:[dev|production]`
Build and export the project as a series of .html files.

### `npm run deploy:[dev|production]`
Export, build and deploy the docker image for release.

## dotenv
We use dotenv / next-dotenv to organize loading of environment variables--more specific settings will override less specific ones. i.e. `.env.development` settings take precedent over `.env`. You should define your own `.env.*.local` which are hidden by git (see [dotenv-load](https://github.com/formatlos/dotenv-load)).

Furthermore, environment prefixes are important for NextJS.

- `NEXT_PUBLIC_*`: available server / client side
- `NEXT_SERVER_*`: available server side only
- `NEXT_STATIC_*`: available during static rendering

Other variables may not propagate, so ensure you use these prefixes.

## UI SCHEMAS
Modifying SigCom UI to your data is done via the entries in the schema table. SigCom uses three types of schemas.

- `/dcic/signature-commons-schema/v5/meta/schema/landing-ui.json`: Defines the overall look of the landing page
- `/dcic/signature-commons-schema/v5/meta/schema/ui-schema.json`: Used for formatting the labels for searches, without this, SigCom defaults to the ids as labels
- `/dcic/signature-commons-schema/v5/meta/schema/counting.json`: Tells sigcom which metadata to count for the landing page, you can also set it to either be a pie chart, a bar chart, or a regular count|

For more information check [Modifying UI](./components/Landing/README.md)
and the [examples folder](./examples/). Also, check [UI_values](./util/ui_values.js)# Modifying the landing page

We can adapt the landing page for signature commons by modifying two json files in the ui-schemas [landing-ui.json](../../examples/dashboard/landing_ui.json) and [ui.json](../../examples/dashboard/ui.json).

![alt text](../../static/sigcom-landing.png)

## Setting fields to count
We can tell the UI which keys to count by creating a database entry on the schema table with this validator: ```/dcic/signature-commons-schema/v5/meta/schema/counting.json```

The following table describe the fields for each entry

| Field         | Value           | Remarks |
| ------------- |---------------| ----------|
| Field_Name     | string          | Name of the meta field in the database|
| Type           | string          | type of the field in the database (string, object)|
| Table | string | Name of the field's table |
| Preferred_Name | string          |Display name of the field|
| Preferred_Name_Singular | string          |Used for pie chart captions|
| Slice | int          |Used to tell how many slices to use for the piecharts|
| MDI_Icon | string | Display icon (see [mdi-icon](https://materialdesignicons.com/) for more information) |
| Meta_Count | boolean | Tells the UI to display the field as part of the meta counts |
| Pie_Count | boolean | Tells the UI to display the field as part of the pie charts |
| Bar_Count | boolean | Tells the UI to display the field as part of the bar charts |
| Table_Count | boolean | Tells the UI to display this as part of the table counts |
| Visible_On_Landing | boolean | If Table_Count is true and this is true, the UI will display the stat on the landing page|
| Visible_On_Admin | boolean | If Table_Count is true and this is true, the UI will display the stat on the Admin page|

Sample entries can be found on [examples folder](../../examples/dashboard). Relevant files are named as ```counting_ui_*.json```. Note that each object in the list is an entry.
```
{
    id: some-uuid-here-to-display,
    meta: {
        # Place values here
    }
}
```

## Changing text content in the page
We can tell the UI which keys to count by creating a database entry on the schema table with this validator: ```/dcic/signature-commons-schema/v5/meta/schema/landing-ui.json```

The following table describe the fields for each entry

| Field         | Value           | Remarks |
| ------------- |---------------| ----------|
| landing | boolean | Use this ui for the landing page|
| admin | boolean | Use this ui for the admin page|
| content | object | You place your modifications here (see below) |

### Content field
Refer to the image above for more information.

| Field         | Value           | Remarks |
| ------------- |---------------| ----------|
| library_name | string | The field that we'll use as the library name (required) |
| resource_name | string | The field that we'll use as the resource name |
| resource_name_from_library | string | Alternatively, if resource table is empty, we can use this field to define which meta field in the library to use. Otherwise, we'll use the library's dataset as resource name |
| resource_icon | string | location of the icon in the metadata |
| signature_search | boolean | tells the UI whether signature_search is activated |
| metadata_search | boolean | tells the UI whether metadata_search is activated |
| resources | boolean | tells the UI whether resources is activated |
| resource_list_style | object | JSS styling for resource page |
| header | string | Header of the UI (Default: Signature Commons) |
| metadata_placeholder | string | Placeholder for the metadata search box |
| geneset_placeholder | string | Placeholder for the geneset search box |
| up_genes_placeholder | string | Placeholder for the up genes search box |
| down_genes_placeholder | string | Placeholder for the down genes search box |
| search_terms | array | Chip terms for metadata search |
| geneset_terms | string | Tab-delimited gene terms for geneset search |
| weighted_geneset_terms | string | Tab-delimited gene weighted terms for geneset search |
| up_set_terms | string | Tab-delimited gene terms for up geneset |
| down_set_terms | string | Tab-delimited gene terms for down geneset search |
| text_1 | string | Text for text_1 |
| text_2 | string | Text for text_2 |
| text_3 | string | Text for text_3 |
| text_4 | string | Text for text_4 |
| resource_pie_caption | string | Caption for resource pie chart |
| bar_chart | object | Controls the barcharts |
| counting_validator | string | Name of the counting validator to use |
| ui_schema| string | name of ui_schema validator to use |
|deactivate_download | boolean | deactivate downloads for this instance|
|deactivate_wordcloud | boolean | deactivate word cloud for this instance|
| bar_chart_style | object | props to pass to rechart (see rechart, ../../util/ui_values.js and ../../examples/dashboard for examples)|
| pie_chart_style | object | props to pass to rechart (see rechart, ../../util/ui_values.js and ../../examples/dashboard for examples)|
|maxResourcesToShow|int|The number of resources triggering More/Less activation at all|
|maxResourcesToShow|number|The maximum resources to show before collapsing if More/Less is activated|
The bar_chart an object with the following fields

| Field         | Value           | Remarks |
| ------------- |---------------| ----------|
| Field_Name | string | Name of the field to create a bar chart (recall that we have a Bar_Count field in the landing-ui.json, we choose which field to construct a bar graph with this field)|
| Caption | string | Caption of the bargraph |

Sample entries can be found on [examples folder](../../examples/dashboard). Relevant files are named as ```ui_*.json```. Note that each object in the list is an entry.
```
{
    id: some-uuid-here-to-display,
    meta: {
        # Place values here
    }
}
```
# Signature Commons UI: Standalone Components

These components are developed for and used by signature commons but are designed to be completely standalone.
# signature-commons-ui-running-sum
A standalone running sum React Component for signature commons
# react-datatable
Datatable components# react-fairshake-insignia
A standalone fairshake insignia react component for sigcom# text-field-suggest
A standalone react component for text-field-sugges

PROPTYPES:
```
TextFieldSuggest.propTypes = {
    input: PropTypes.arrayOf(PropTypes.shape({
        label: PropTypes.string.isRequired,
        type: PropTypes.oneOf(["valid", "suggestions", "invalid", "loading", "disabled"]),
        id: PropTypes.oneOfType([
            PropTypes.string,
            PropTypes.number
          ]),
        suggestions: PropTypes.arrayOf(PropTypes.shape({
            label: PropTypes.string.isRequired,
            type: PropTypes.oneOf(["valid", "suggestions", "invalid", "loading", "disabled"]),
            id: PropTypes.oneOfType([
                PropTypes.string,
                PropTypes.number
              ]),     
        })),
        gridColumnProps: PropTypes.object,
        gridRowProps: PropTypes.object,
        avatarProps: PropTypes.object,
        labelProps: PropTypes.object,
        chipProps: PropTypes.object,
        suggestionsProps: PropTypes.object,
    })).isRequired,
    onSubmit: PropTypes.func.isRequired,
    onAdd: PropTypes.func.isRequired,
    onDelete: PropTypes.func.isRequired,
    onClick: PropTypes.func.isRequired,
    onSuggestionClick: PropTypes.func.isRequired,
    renderChip: PropTypes.func,
    colors_and_icon: PropTypes.shape({
        background: PropTypes.string,
        color: PropTypes.string,
        icon: PropTypes.string
    }),
    placeholder: PropTypes.string,
    gridColumnProps: PropTypes.object,
    gridRowProps: PropTypes.object,
    avatarProps: PropTypes.object,
    labelProps: PropTypes.object,
    chipProps: PropTypes.object,
    chipInputProps: PropTypes.object,
    formProps: PropTypes.object,
    suggestionsProps: PropTypes.object,
}
``````
{
  fetch_meta: true
  count: true,
  per_parent_count: true,
  value_count: true
  query: {
    search: [...],
    filters: {
      [filter_field]: [...]
    },
    skip,
    limit
  },
  parent_ids: [...]
  value_count_params: {
    fields: [...],
    skip,
    limit
  }
}
```

```
url params
{
  search: [...],
  libraries: {
    filters: {
      [filter_field]: [...]
    },
    skip,
    limit,
    order,
    value_count_params: {
        fields: [...]
    }
  },
  signatures: {
    filters: {
      [filter_field]: [...]
    },
    skip,
    limit
  }  
}
```

```
saga params
{
  search: [...],
  libraries: {
    filters: {
      [filter_field]: [...]
    },
    skip,
    limit,
    operations: {
      metadata_search: false,
      per_parent_count: false,
      value_count: false,
      count: true
    }
  },
  signatures: {
    filters: {
      [filter_field]: [...]
    },
    skip,
    limit,
    order,
    operations: {
      metadata_search: true,
      per_parent_count: true,
      value_count: true,
      count: true
    }  
  }
}
```# UI Schemas

This will be merged into `signature-commons-schema` when it becomes possible, for now it needs to exist to design the UI around the idea that they would have been present in the schema.
