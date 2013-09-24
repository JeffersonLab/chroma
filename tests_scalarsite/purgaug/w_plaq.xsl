<?xml version="1.0"?>
<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>

<xsl:output method="text"/>
<xsl:template match="/">

<xsl:for-each select="//purgaug/doHB/MCUpdates/elem/Update/InlineObservables/elem/Plaquette">

  <xsl:value-of select="w_plaq" />
  <xsl:text>
  </xsl:text>
</xsl:for-each>

</xsl:template>


</xsl:stylesheet> 
