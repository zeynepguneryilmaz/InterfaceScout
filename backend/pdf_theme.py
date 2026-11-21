"""Consistent, restrained print styling for InterfaceScout research PDFs.

Uses real viewer images and existing analysis values; never fabricates visuals,
numerical color legends or docking/adsorption interpretations.
"""
from __future__ import annotations

from io import BytesIO
from colorsys import hls_to_rgb
from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import Paragraph, Image, Spacer, Table, TableStyle, SimpleDocTemplate
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.lib.utils import ImageReader
from reportlab.graphics.shapes import Drawing, Rect, String, Line

BURGUNDY = colors.HexColor("#6B4355")
SAGE = colors.HexColor("#718375")
INK = colors.HexColor("#302C34")
MUTED = colors.HexColor("#6F6871")
LINE = colors.HexColor("#E4DFE1")
IVORY = colors.HexColor("#FAF8F5")
PAGE_W, PAGE_H = A4

_styles = getSampleStyleSheet()
_styles.add(ParagraphStyle("AtlasEyebrow", parent=_styles["Normal"], fontName="Helvetica-Bold",
    fontSize=8, leading=11, textColor=BURGUNDY, spaceAfter=5))
_styles.add(ParagraphStyle("AtlasTitle", parent=_styles["Normal"], fontName="Helvetica-Bold",
    fontSize=21, leading=25, textColor=INK, spaceAfter=9))
_styles.add(ParagraphStyle("AtlasHeading", parent=_styles["Normal"], fontName="Helvetica-Bold",
    fontSize=12, leading=16, textColor=BURGUNDY, spaceBefore=11, spaceAfter=7))
_styles.add(ParagraphStyle("AtlasBody", parent=_styles["Normal"], fontName="Helvetica",
    fontSize=9, leading=13, textColor=INK, spaceAfter=7, allowWidows=0, allowOrphans=0))
_styles.add(ParagraphStyle("AtlasCaption", parent=_styles["Normal"], fontName="Helvetica",
    fontSize=8, leading=11, textColor=MUTED, spaceAfter=7))
_styles.add(ParagraphStyle("AtlasCell", parent=_styles["Normal"], fontName="Helvetica",
    fontSize=8, leading=10, textColor=INK))
_styles.add(ParagraphStyle("AtlasCellHead", parent=_styles["AtlasCell"], fontName="Helvetica-Bold", textColor=BURGUNDY))

def paragraph(message, kind="AtlasBody"):
    return Paragraph(escape(str(message if message is not None else "")), _styles[kind])

def heading(text):
    return paragraph(text,"AtlasHeading")

def table(headers, rows, widths=None):
    data = [[paragraph(h,"AtlasCellHead") for h in headers]]
    data.extend([[paragraph(value,"AtlasCell") for value in row] for row in rows])
    if len(data)==1:
        data.append([paragraph("No ranked regions","AtlasCell")]+[paragraph("","AtlasCell") for _ in headers[1:]])
    t=Table(data,colWidths=widths,repeatRows=1,hAlign="LEFT",splitByRow=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),IVORY),
        ("LINEBELOW",(0,0),(-1,0),0.8,LINE),
        ("LINEBELOW",(0,1),(-1,-1),0.25,LINE),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("LEFTPADDING",(0,0),(-1,-1),8),
        ("RIGHTPADDING",(0,0),(-1,-1),8),
        ("TOPPADDING",(0,0),(-1,-1),7),
        ("BOTTOMPADDING",(0,0),(-1,-1),7),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,IVORY]),
    ]))
    return t

def viewer_figure(blob: bytes, max_width=495, max_height=292):
    im=PILImage.open(BytesIO(blob))
    if im.width<=0 or im.height<=0:
        raise ValueError("Invalid viewer screenshot")
    scale=min(max_width/im.width,max_height/im.height)
    return Image(BytesIO(blob),width=im.width*scale,height=im.height*scale)

def page_decoration(canvas,doc):
    canvas.saveState()
    canvas.setFillColor(IVORY)
    canvas.rect(0,PAGE_H-39,PAGE_W,39,fill=1,stroke=0)
    logo=Path(__file__).resolve().parents[1]/"frontend"/"logo.png"
    if logo.exists():
        try:
            canvas.drawImage(ImageReader(str(logo)),37,PAGE_H-34,26,26,mask="auto")
        except Exception:
            pass
    canvas.setFillColor(INK)
    canvas.setFont("Helvetica-Bold",10)
    canvas.drawString(71,PAGE_H-25,"InterfaceScout")
    canvas.setFillColor(BURGUNDY)
    canvas.setFont("Helvetica-Bold",7.5)
    canvas.drawRightString(PAGE_W-37,PAGE_H-24,"SURFACE ATLAS")
    canvas.setStrokeColor(LINE)
    canvas.line(37,40,PAGE_W-37,40)
    canvas.setFillColor(MUTED)
    canvas.setFont("Helvetica",7.2)
    canvas.drawString(37,27,"Descriptive physicochemical mapping · Not an adsorption prediction")
    canvas.drawRightString(PAGE_W-37,27,str(doc.page))
    canvas.restoreState()


def consensus_frequency_key():
    """Exact HSL transformation used by frontend. Frequency 0 is not rendered."""
    values=[i/10 for i in range(1,11)]
    swatches=[]
    for f in values:
        h=(228-207*f)/360
        light=.49+.08*f
        red,green,blue=hls_to_rgb(h,light,.82)
        swatches.append(colors.Color(red,green,blue))
    labels=["0.1", "", "", "", "0.5", "", "", "", "", "1.0"]
    tab=Table([[""]*10,labels],colWidths=[50.4]*10,rowHeights=[12,13],hAlign="LEFT")
    style=[
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
        ("ALIGN",(0,0),(-1,-1),"CENTER"),
        ("TEXTCOLOR",(0,1),(-1,1),MUTED),
        ("FONTSIZE",(0,1),(-1,1),7),
    ]
    for i,c in enumerate(swatches):
        style.append(("BACKGROUND",(i,0),(i,0),c))
    tab.setStyle(TableStyle(style))
    return tab


def jaccard_chart(map_comparisons):
    """Vector PDF chart of Top-5 Jaccard, a bounded descriptive overlap score."""
    items=[]
    for m in map_comparisons:
        ks=m.get("top_k") or []
        if len(ks)<3:
            continue
        value=max(0,min(1,float(ks[2].get("jaccard") or 0)))
        items.append((str(m.get("key","")),value))
    width=504
    row_h=20
    height=40+row_h*len(items)
    drawing=Drawing(width,height)
    bar_x=160; bar_w=292
    for k,tick in enumerate((0,.5,1.0)):
        x=bar_x+tick*bar_w
        drawing.add(Line(x,25,x,height-15,strokeColor=LINE,strokeWidth=.8))
        drawing.add(String(x,12,str(tick),fontName="Helvetica",fontSize=7,
                           fillColor=MUTED,textAnchor="middle"))
    for i,(name,value) in enumerate(items):
        y=height-33-i*row_h
        drawing.add(String(1,y+3,name[:23],fontName="Helvetica",fontSize=7.5,
                           fillColor=INK))
        drawing.add(Rect(bar_x,y-1,bar_w,12,fillColor=IVORY,strokeColor=None))
        if value>0:
            drawing.add(Rect(bar_x,y-1,bar_w*value,12,fillColor=SAGE,strokeColor=None))
        drawing.add(String(bar_x+bar_w+8,y+2,f"{value:.2f}",
                           fontName="Helvetica-Bold",fontSize=8,fillColor=BURGUNDY))
    return drawing

class AtlasDocTemplate(SimpleDocTemplate):
    """Page references and PDF outline entries from real section flowables."""
    def afterFlowable(self, flowable):
        key=getattr(flowable,"atlas_bookmark",None)
        if not key:return
        level=getattr(flowable,"atlas_level",0)
        title=getattr(flowable,"atlas_title","Section")
        self.canv.bookmarkPage(key)
        self.canv.addOutlineEntry(title,key,level=level,closed=False)
        if getattr(flowable,"atlas_in_toc",True):
            self.notify("TOCEntry",(level,title,self.page,key))


def atlas_section(text,key,level=0,toc=True):
    obj=paragraph(text,"AtlasHeading" if level else "AtlasTitle")
    obj.atlas_bookmark=key
    obj.atlas_title=text
    obj.atlas_level=level
    obj.atlas_in_toc=toc
    return obj


def atlas_contents():
    toc=TableOfContents()
    toc.levelStyles=[
        ParagraphStyle("AtlasTOC0",fontName="Helvetica",fontSize=10,leading=17,
                       leftIndent=6,rightIndent=16,spaceBefore=3,textColor=INK),
        ParagraphStyle("AtlasTOC1",fontName="Helvetica",fontSize=9,leading=15,
                       leftIndent=18,rightIndent=16,spaceBefore=2,textColor=MUTED),
    ]
    return toc
