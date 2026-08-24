import json
import base64
import io
import textwrap

from PIL import Image as PILImage
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, HRFlowable, Preformatted
)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

TTF_MONO = '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'
pdfmetrics.registerFont(TTFont('DejaVuMono', TTF_MONO))

NOTEBOOK_PATH = '1bis_Molecule_Analysis.ipynb'

# ── load notebook ──────────────────────────────────────────────────────────────
with open(NOTEBOOK_PATH) as f:
    nb = json.load(f)
cell_map = {c.get('id', ''): c for c in nb['cells']}

# ── helpers ────────────────────────────────────────────────────────────────────
def get_identifier():
    src = ''.join(cell_map['config-cell']['source'])
    for line in src.split('\n'):
        if line.startswith('IDENTIFIER'):
            return line.split('=')[1].strip().strip("'\"")
    return 'unknown'

def cell_png_flowable(cell, max_width=14*cm, max_height=10*cm):
    for out in cell.get('outputs', []):
        if out.get('output_type') in ('display_data', 'execute_result'):
            png = out.get('data', {}).get('image/png', '')
            if png:
                if isinstance(png, list):
                    png = ''.join(png)
                img_bytes = base64.b64decode(png)
                pil = PILImage.open(io.BytesIO(img_bytes))
                w, h = pil.size
                ratio = min(max_width / w, max_height / h)
                return Image(io.BytesIO(img_bytes), width=w*ratio, height=h*ratio)
    return None

def cell_text(cell):
    parts = []
    for out in cell.get('outputs', []):
        if out.get('output_type') == 'stream':
            t = out.get('text', '')
            parts.append(''.join(t) if isinstance(t, list) else t)
        elif out.get('output_type') in ('display_data', 'execute_result'):
            t = out.get('data', {}).get('text/plain', '')
            if isinstance(t, list):
                t = ''.join(t)
            if t and not t.startswith('<PIL') and not t.startswith('<Figure'):
                parts.append(t)
    return ''.join(parts).rstrip()

# ── styles ─────────────────────────────────────────────────────────────────────
styles = getSampleStyleSheet()
style_title = ParagraphStyle('Title',
    fontSize=28, leading=36, spaceAfter=6,
    textColor=colors.HexColor('#1a1a2e'), alignment=1, fontName='Helvetica-Bold')
style_subtitle = ParagraphStyle('Subtitle',
    fontSize=14, leading=20, spaceAfter=20,
    textColor=colors.HexColor('#555555'), alignment=1)
style_h2 = ParagraphStyle('H2',
    fontSize=13, leading=18, spaceBefore=14, spaceAfter=6,
    textColor=colors.HexColor('#1a1a2e'), fontName='Helvetica-Bold')
style_mono = ParagraphStyle('Mono',
    fontSize=9, leading=13, fontName='DejaVuMono',
    textColor=colors.HexColor('#222222'))
style_mono_small = ParagraphStyle('MonoSmall',
    fontSize=8, leading=12, fontName='DejaVuMono',
    textColor=colors.HexColor('#333333'))
hr = lambda: HRFlowable(width='100%', thickness=0.5,
                         color=colors.HexColor('#cccccc'), spaceAfter=8)

# ── build document ─────────────────────────────────────────────────────────────
identifier = get_identifier()
OUTPUT_PDF = f'{identifier}.pdf'

doc = SimpleDocTemplate(OUTPUT_PDF, pagesize=A4,
                        leftMargin=2*cm, rightMargin=2*cm,
                        topMargin=2.5*cm, bottomMargin=2*cm)

story = []

# title
story.append(Spacer(1, 3*cm))
story.append(Paragraph(identifier, style_title))
story.append(Paragraph('Single Molecule Analysis', style_subtitle))
story.append(hr())

# 2D structure
story.append(Paragraph('2D Structure', style_h2))
img_flowable = cell_png_flowable(cell_map['structure-cell'])
if img_flowable:
    story.append(img_flowable)
    story.append(Spacer(1, 0.3*cm))
smiles_inchi = cell_text(cell_map['structure-cell'])
if smiles_inchi:
    wrapped = '\n'.join(
        textwrap.fill(line, width=95) if len(line) > 95 else line
        for line in smiles_inchi.splitlines()
    )
    story.append(Preformatted(wrapped, style_mono_small))

story.append(Spacer(1, 0.4*cm))
story.append(hr())

# molecular properties
story.append(Paragraph('Molecular Properties', style_h2))
story.append(Preformatted(cell_text(cell_map['mol-cell']), style_mono))

story.append(Spacer(1, 0.4*cm))
story.append(hr())

# crystallographic info
story.append(Paragraph('Crystallographic Info', style_h2))
story.append(Preformatted(cell_text(cell_map['crystallo-cell']), style_mono))

doc.build(story)
print(f'PDF generated → {OUTPUT_PDF}')
