<?xml version="1.0"?>
<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>
<xsl:output method="html"/>


<xsl:template match="/testResults">
<html>
<head>
</head>
<body>
<h2>Regression test results</h2>
<p><em> Number of tests: <xsl:value-of select="//summary/@tests"/></em></p>
<p><em> Number of successes: <xsl:value-of select="//summary/@successes"/></em></p>
<p><em> Number of failures:<xsl:value-of select="//summary/@failures"/></em></p> 
<table border='2'>		
<tr>	
<td><b>PROGRAM</b></td> <td><b>Candidate</b></td> <td><b>Result</b></td>
</tr>

<xsl:for-each select="./test">
    <xsl:choose>
	<xsl:when test="string(.)='PASS'">
		<tr bgcolor="#99ff33">
		<td>		
			<xsl:value-of select="@program"/>
		</td>
		<td>
			<xsl:value-of select="@candidate"/>
		</td>
		<td>
			<xsl:value-of select="." />
		</td>
		</tr>
	</xsl:when>
    	<xsl:otherwise>
		<tr bgcolor="#ff3300">
		<td text="#ffff33">		
			<xsl:value-of select="@program"/>
		</td>
		<td>
			<xsl:value-of select="@candidate"/>
		</td>
		<td>
			<xsl:value-of select="." />
		</td>
		</tr>
	</xsl:otherwise>
	</xsl:choose>	

</xsl:for-each>
</table>
</body>
</html>
</xsl:template>


</xsl:stylesheet>
