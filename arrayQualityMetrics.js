// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, true, true, true, true, true, true ];
var arrayMetadata    = [ [ "1", "GSM800742_chip_array_C06N-H.CEL.gz", "1", "08/01/08 12:23:56" ], [ "2", "GSM800743_chip_array_C11N-H.CEL.gz", "2", "08/01/08 12:35:17" ], [ "3", "GSM800744_chip_array_C24N-H.CEL.gz", "3", "08/07/08 17:15:10" ], [ "4", "GSM800745_chip_array_C27N-H.CEL.gz", "4", "08/01/08 13:52:28" ], [ "5", "GSM800746_chip_array_C28N-H.CEL.gz", "5", "08/01/08 14:03:46" ], [ "6", "GSM800747_chip_array_C30N-H.CEL.gz", "6", "08/01/08 15:00:05" ], [ "7", "GSM800748_chip_array_C31N-H.CEL.gz", "7", "08/01/08 15:11:16" ], [ "8", "GSM800749_chip_array_C32N-H.CEL.gz", "8", "08/01/08 15:22:42" ], [ "9", "GSM800750_chip_array_C33N-H.CEL.gz", "9", "08/01/08 15:33:53" ], [ "10", "GSM800751_chip_array_C35N-H.CEL.gz", "10", "08/07/08 12:17:15" ], [ "11", "GSM800752_chip_array_C36N-H.CEL.gz", "11", "08/07/08 12:29:11" ], [ "12", "GSM800753_chip_array_C38N-H.CEL.gz", "12", "08/07/08 12:40:31" ], [ "13", "GSM800754_chip_array_C41N-H.CEL.gz", "13", "08/07/08 12:51:39" ], [ "14", "GSM800755_chip_array_C42N-H.CEL.gz", "14", "08/07/08 13:48:28" ], [ "15", "GSM800756_chip_array_C44N-H.CEL.gz", "15", "08/07/08 13:59:57" ], [ "16", "GSM800757_chip_array_C45N-H.CEL.gz", "16", "08/07/08 14:11:08" ], [ "17", "GSM800758_chip_array_C47N-H.CEL.gz", "17", "07/12/07 12:09:45" ], [ "18", "GSM800759_chip_array_C06T-H.CEL.gz", "18", "10/18/07 14:39:39" ], [ "19", "GSM800760_chip_array_C11T-H.CEL.gz", "19", "10/18/07 14:50:56" ], [ "20", "GSM800761_chip_array_C24T-H.CEL.gz", "20", "10/18/07 16:09:16" ], [ "21", "GSM800762_chip_array_C27T-H.CEL.gz", "21", "10/18/07 16:32:30" ], [ "22", "GSM800763_chip_array_C28T-H.CEL.gz", "22", "10/18/07 16:43:48" ], [ "23", "GSM800764_chip_array_C30T-H.CEL.gz", "23", "10/19/07 12:00:32" ], [ "24", "GSM800765_chip_array_C31T-H.CEL.gz", "24", "10/19/07 12:11:51" ], [ "25", "GSM800766_chip_array_C32T-H.CEL.gz", "25", "10/19/07 12:22:53" ], [ "26", "GSM800767_chip_array_C33T-H.CEL.gz", "26", "10/19/07 12:34:04" ], [ "27", "GSM800768_chip_array_C35T-H.CEL.gz", "27", "10/19/07 13:26:46" ], [ "28", "GSM800769_chip_array_C36T-H.CEL.gz", "28", "10/19/07 13:38:32" ], [ "29", "GSM800770_chip_array_C38T-H.CEL.gz", "29", "10/19/07 13:50:06" ], [ "30", "GSM800771_chip_array_C41T-H.CEL.gz", "30", "10/19/07 14:01:29" ], [ "31", "GSM800772_chip_array_C42T-H.CEL.gz", "31", "10/19/07 14:57:14" ], [ "32", "GSM800773_chip_array_C44T-H.CEL.gz", "32", "10/19/07 15:08:27" ], [ "33", "GSM800774_chip_array_C45T-H.CEL.gz", "33", "10/19/07 15:19:46" ], [ "34", "GSM800775_chip_array_C47T-H.CEL.gz", "34", "10/19/07 15:31:03" ], [ "35", "GSM800776_chip_array_colo320_0.5AZA.CEL.gz", "35", "12/09/09 13:56:16" ], [ "36", "GSM800777_chip_array_colo320_control.CEL.gz", "36", "12/09/09 13:45:00" ], [ "37", "GSM800778_chip_array_HCT116_0.5AZA.CEL.gz", "37", "12/09/09 11:15:47" ], [ "38", "GSM800779_chip_array_HCT116_control.CEL.gz", "38", "12/09/09 11:04:22" ], [ "39", "GSM800780_chip_array_HT29_0.5AZA.CEL.gz", "39", "12/09/09 12:18:06" ], [ "40", "GSM800781_chip_array_HT29_control.CEL.gz", "40", "12/09/09 12:06:21" ], [ "41", "GSM800782_chip_array_RKO_0.5AZA.CEL.gz", "41", "12/09/09 10:52:55" ], [ "42", "GSM800783_chip_array_RKO_control.CEL.gz", "42", "12/09/09 10:41:22" ], [ "43", "GSM800784_chip_array_SW480_0.5AZA.CEL.gz", "43", "12/09/09 12:40:30" ], [ "44", "GSM800785_chip_array_SW480_control.CEL.gz", "44", "12/09/09 12:29:23" ] ];
var svgObjectNames   = [ "pca", "dens", "dig" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
