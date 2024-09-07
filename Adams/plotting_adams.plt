$---------------------------------------------------------------------MDI_HEADER
[MDI_HEADER]
 FILE_TYPE     =  'plt'
 FILE_VERSION  =  2.0
 FILE_FORMAT   =  'ASCII'
$---------------------------------------------------------------------------PAGE
[PAGE]
 PAGE_LAYOUT  =  11.0
 NUMBER_OF_PLOTS  =  1.0
 PAGE_NAME  =  'page_1'
$---------------------------------------------------------------------------PLOT
[PLOT]
 INDEX  =  0.0
 NAME  =  'plot_1'
 TIME_LOWER_LIMIT  =  0.0
 TIME_UPPER_LIMIT  =  0.0
(LEGEND)
{placement   location     fill   grow_left   grow_down   font}
'top left'   2.1,81.6     1      FALSE       TRUE        9
(PLOT_BORDER)
{color   line_style    line_weight}
'BLACK'  'solid'       1.0
(PRIMARY_GRID)
{color    line_style    line_weight}
'SILVER'  'solid'       0.5
(SECONDARY_GRID)
{color    line_style    line_weight}
'SILVER'  'solid'       0.5
(LEGEND_BORDER)
{color   line_style    line_weight}
'BLACK'  'solid'       1.0
(GRAPH_AREA)
{minX        minY        maxX        maxY        auto_graph_area}
 -0.7572      12.7147   160.7572     83.7419     yes
(NOTES)
{name       type             color     placement      alignment         location     font     autopos     autogenerate     numStrings}
'header'    'table header'   'BLACK'   'horizontal'   'center_bottom'   8.4,1.9      1        no          yes              2
 STRING_1_TEXT  =  'adams_project'
 STRING_2_TEXT  =  'JOINT_C: JOINT/4'
{name         type         color     placement      alignment      location     font     autopos     autogenerate     numStrings}
'analysis'    'analysis'   'BLACK'   'horizontal'   'center_top'   -0.8,5.5     9        no          yes              1
 STRING_1_TEXT  =  'Analysis:  Last_Run'
{name     type     color     placement      alignment      location      font     autopos     autogenerate     numStrings}
'date'    'date'   'BLACK'   'horizontal'   'center_top'   160.8,5.5     9        no          yes              1
 STRING_1_TEXT  =  '2021-07-16 10:02:10'
{name      type      color     placement      alignment         location      font     autopos     autogenerate     numStrings}
'title'    'title'   'BLACK'   'horizontal'   'center_bottom'   80.0,88.8     9        no          yes              1
 STRING_1_TEXT  =  'adams_project'
{name         type         color     placement      alignment         location      font     autopos     autogenerate     numStrings}
'subtitle'    'subtitle'   'BLACK'   'horizontal'   'center_bottom'   80.0,85.3     9        no          yes              1
 STRING_1_TEXT  =  'JOINT_C: JOINT/4'
(PLOT_AXES_FORMAT)
{axis_name    type          color    placement    scaling    offset    primary    limits}
'vaxis'       'vertical'    'BLACK'  'left'       'linear'   0.0       yes        0.000000,0.000000
'haxis'       'horizontal'  'BLACK'  'bottom'     'linear'   0.0       yes        0.000000,0.000000
(PLOT_AXES_LABELS)
{axis_name    label             color     placement     alignment        font    autopos    offset    location}
'vaxis'       'Force (newton)'  'BLACK'   'vertical'    'center_bottom'  9       1          10.7      -11.5,48.2
'haxis'       'Time (sec)'      'BLACK'   'horizontal'  'center_top'     9       1          7.2       80.0,5.5
(PLOT_AXES_TICS)
{axis_name    auto_divisions    use_divisions    divisions    increments    minor_divisions    color}
'vaxis'       'yes'             'yes'            4            5.000         2                  'BLACK'
'haxis'       'yes'             'yes'            6            0.500         2                  'BLACK'
(PLOT_AXES_NUMBERS)
{axis_name    trailing_zeros    decimal_places    scientific_range    font    color}
'vaxis'       0                 4                 -4,5                9.0     'BLACK'
'haxis'       0                 4                 -4,5                9.0     'BLACK'
$---------------------------------------------------------------------PLOT_CURVE
[PLOT_CURVE]
 NAME  =  'curve_1'
 PLOT  =  'plot_1'
 VERTICAL_AXIS  =  'vaxis'
 HORIZONTAL_AXIS  =  'haxis'
 HORIZONTAL_EXPRESSION  =  'JOINT_C.TIME'
 HORIZONTAL_COMPONENT  =  'JOINT_C.TIME'
 VERTICAL_EXPRESSION  =  'JOINT_C.FX'
 VERTICAL_COMPONENT  =  'JOINT_C.FX'
 Y_UNITS  =  'force'
 X_UNITS  =  'time'
 LEGEND_TEXT  =  'JOINT_C.FX'
 COLOR  =  'red'
 STYLE  =  'solid'
 SYMBOL  =  'NONE'
 LINE_WEIGHT  =  2.0
 HOTPOINT  =  0.0
 INCREMENT_SYMBOL  =  1.0
