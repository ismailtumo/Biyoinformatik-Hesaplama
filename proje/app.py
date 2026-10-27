import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Restriction import EcoRI, BamHI, HindIII, NotI, XhoI, TaqI
from Bio.Blast import NCBIWWW, NCBIXML
import io
import random

# Sayfa Ayarları
st.set_page_config(page_title="Gelişmiş Biyoinformatik Analizörü", layout="wide")

# --- STREAMLIT ARAYÜZ ELEMANLARINI GİZLEME (GÜÇLENDİRİLMİŞ CSS) ---
hide_streamlit_style = """
<style>
/* Sağ üstteki Fork, GitHub butonu ve üç noktayı gizler */
[data-testid="stToolbar"] {visibility: hidden !important; display: none !important;}
[data-testid="stHeader"] {visibility: hidden !important; display: none !important;}
header {visibility: hidden !important; display: none !important;}

/* Alt kısımdaki yazıları ve logoları gizler */
footer {visibility: hidden !important;}
#MainMenu {visibility: hidden !important;}

/* Sağ alttaki geliştirici menüsü / rozetleri gizler */
.viewerBadge_container {display: none !important;}
.stDeployButton {display: none !important;}
</style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.title("🧬 Kapsamlı DNA Analiz ve Gen Anotasyon Aracı")
st.markdown("""
Bu araç, DNA dizilerinin biyoinformatik analizlerini yapar, Açık Okuma Çerçevelerini (ORF) bulur, 
protein biyokimyasını hesaplar, enzim kesim bölgelerini haritalar ve NCBI BLAST ile anotasyon sağlar.
""")

if 'ornek_dizi' not in st.session_state:
    st.session_state['ornek_dizi'] = ""

st.sidebar.header("1. Veri Yükleme")

if st.sidebar.button("🧪 Rastgele Örnek Dize Üret (1.5kb)"):
    st.session_state['ornek_dizi'] = "".join(random.choices(["A", "C", "G", "T"], k=1500))

input_method = st.sidebar.radio("Giriş Yöntemi Seçin:", ("Metin Girişi", "FASTA/GenBank Dosyası"))

dna_sequence = ""

if input_method == "Metin Girişi":
    raw_seq = st.sidebar.text_area("DNA Dizisini Buraya Yapıştırın:", value=st.session_state['ornek_dizi'], height=150)
    dna_sequence = "".join(raw_seq.split()).upper()
else:
    uploaded_file = st.sidebar.file_uploader("Bir dizi dosyası seçin", type=["fasta", "txt", "fa", "gb"])
    if uploaded_file is not None:
        stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        lines = stringio.readlines()
        # Basit FASTA okuyucu
        dna_sequence = "".join([line.strip().upper() for line in lines if not line.startswith(">")])

# --- FONKSİYONLAR ---
def validate_dna(seq):
    valid_chars = set("ACGTN")
    return all(char in valid_chars for char in set(seq))

def calculate_protein_properties(protein_seq):
    """Proteinin molekül ağırlığı ve izoelektrik noktasını hesaplar"""
    clean_seq = protein_seq.replace("*", "") # Durdurma kodonunu çıkar
    if "X" in clean_seq or "B" in clean_seq or "Z" in clean_seq:
        return None, None # Belirsiz amino asitler varsa hesaplama yapılamaz
    
    try:
        analysis = ProteinAnalysis(clean_seq)
        mw = analysis.molecular_weight() / 1000 # kDa cinsinden
        pi = analysis.isoelectric_point()
        return round(mw, 2), round(pi, 2)
    except:
        return None, None

def find_orfs(seq, min_protein_length=10):
    orfs = []
    seq_len = len(seq)
    
    # İleri çerçeveler (+1, +2, +3)
    for frame in range(3):
        for i in range(frame, seq_len - 2, 3):
            codon = seq[i:i+3]
            if codon == "ATG": 
                for j in range(i+3, seq_len - 2, 3):
                    if seq[j:j+3] in ["TAA", "TAG", "TGA"]: 
                        orf_seq = seq[i:j+3]
                        protein = str(Seq(orf_seq).translate(to_stop=True))
                        if len(protein) >= min_protein_length:
                            mw, pi = calculate_protein_properties(protein)
                            orfs.append({
                                "Yön": "+",
                                "Çerçeve": f"+{frame+1}",
                                "Başlangıç": i + 1,
                                "Bitiş": j + 3,
                                "DNA Uzunluk": len(orf_seq),
                                "AA Uzunluk": len(protein),
                                "Mol. Ağırlık (kDa)": mw,
                                "pI": pi,
                                "Protein Dizisi": protein
                            })
                        break
                        
    # Geri çerçeveler (-1, -2, -3)
    rev_seq = str(Seq(seq).reverse_complement())
    for frame in range(3):
        for i in range(frame, seq_len - 2, 3):
            if rev_seq[i:i+3] == "ATG":
                for j in range(i+3, seq_len - 2, 3):
                    if rev_seq[j:j+3] in ["TAA", "TAG", "TGA"]:
                        orf_seq = rev_seq[i:j+3]
                        protein = str(Seq(orf_seq).translate(to_stop=True))
                        if len(protein) >= min_protein_length:
                            mw, pi = calculate_protein_properties(protein)
                            start_pos = seq_len - (j + 3) + 1
                            end_pos = seq_len - i
                            orfs.append({
                                "Yön": "-",
                                "Çerçeve": f"-{frame+1}",
                                "Başlangıç": start_pos,
                                "Bitiş": end_pos,
                                "DNA Uzunluk": len(orf_seq),
                                "AA Uzunluk": len(protein),
                                "Mol. Ağırlık (kDa)": mw,
                                "pI": pi,
                                "Protein Dizisi": protein
                            })
                        break
    return orfs

def plot_gc_content(seq, window_size=50):
    """Genom boyunca GC içeriğini hesaplar ve çizer"""
    if len(seq) < window_size:
        return None
    
    gc_values = []
    positions = []
    for i in range(0, len(seq) - window_size + 1, 10):
        window = seq[i:i+window_size]
        gc = gc_fraction(window) * 100
        gc_values.append(gc)
        positions.append(i + window_size//2)
        
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(positions, gc_values, color='#2ca02c')
    ax.set_title(f"GC İçeriği Dağılımı (Pencere Boyutu: {window_size} bp)")
    ax.set_xlabel("Nükleotid Pozisyonu")
    ax.set_ylabel("GC %")
    ax.grid(True, linestyle='--', alpha=0.7)
    return fig

def plot_orf_map(orfs, seq_len):
    """ORF'lerin dizi üzerindeki konumunu grafiksel olarak haritalar"""
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.plot([0, seq_len], [0, 0], color='black', linewidth=2) # Ana DNA çizgisi
    
    y_offset = 0.5
    for idx, orf in enumerate(orfs):
        start, end = orf['Başlangıç'], orf['Bitiş']
        direction = orf['Yön']
        color = 'blue' if direction == '+' else 'red'
        y_pos = y_offset if direction == '+' else -y_offset
        
        # Ok çizimi
        dx = end - start if direction == '+' else start - end
        x_start = start if direction == '+' else end
        
        ax.arrow(x_start, y_pos, dx, 0, width=0.1, head_width=0.3, head_length=min(20, abs(dx)), 
                 length_includes_head=True, color=color, alpha=0.7)
        ax.text((start+end)/2, y_pos + 0.2 if direction == '+' else y_pos - 0.3, 
                f"ORF {idx+1}", ha='center', fontsize=8)

    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([])
    ax.set_title("ORF Genom Haritası (Mavi: İleri, Kırmızı: Geri)")
    ax.set_xlabel("Nükleotid Pozisyonu")
    return fig

# --- ANALİZ VE ARAYÜZ ---
if dna_sequence:
    if validate_dna(dna_sequence):
        st.success("DNA Dizisi geçerli formatta!")
        
        # 1. TEMEL İSTATİSTİKLER
        st.header("1. Temel İstatistikler ve Nükleotid Kompozisyonu")
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Toplam Uzunluk", f"{len(dna_sequence)} bp")
        col2.metric("GC İçeriği", f"{gc_fraction(dna_sequence) * 100:.2f} %")
        col3.metric("A/T Sayısı", f'{dna_sequence.count("A")} / {dna_sequence.count("T")}')
        col4.metric("G/C Sayısı", f'{dna_sequence.count("G")} / {dna_sequence.count("C")}')

        # Nükleotid ve GC Grafikleri
        col_plot1, col_plot2 = st.columns(2)
        with col_plot1:
            counts = {"A": dna_sequence.count("A"), "T": dna_sequence.count("T"), 
                      "G": dna_sequence.count("G"), "C": dna_sequence.count("C")}
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.bar(counts.keys(), counts.values(), color=['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e'])
            ax.set_title("Nükleotid Dağılımı")
            st.pyplot(fig)
            
        with col_plot2:
            gc_fig = plot_gc_content(dna_sequence)
            if gc_fig:
                st.pyplot(gc_fig)

        # 2. RESTRİKSİYON ENZİM ANALİZİ
        st.header("2. Restriksiyon Enzimi Kesim Haritası")
        my_seq = Seq(dna_sequence)
        enzymes = [EcoRI, BamHI, HindIII, NotI, XhoI, TaqI]
        cut_sites = []
        for enz in enzymes:
            sites = enz.search(my_seq)
            if sites:
                cut_sites.append({"Enzim": enz.__name__, "Kesim Noktaları": str(sites), "Kesim Sayısı": len(sites)})
        
        if cut_sites:
            st.table(pd.DataFrame(cut_sites))
        else:
            st.info("Seçili yaygın enzimler için bu dizide kesim bölgesi bulunamadı.")

        # 3. ORF ANALİZİ VE HARİTALAMA
        st.header("3. Açık Okuma Çerçeveleri (ORF) Analizi")
        min_len = st.number_input("Minimum Protein Uzunluğu (Amino Asit):", min_value=10, max_value=2000, value=20, step=1)
        
        with st.spinner("ORF'ler aranıyor ve protein özellikleri hesaplanıyor..."):
            orf_results = find_orfs(dna_sequence, min_len)
            
        if orf_results:
            st.success(f"Toplam {len(orf_results)} adet ORF bulundu.")
            df_orfs = pd.DataFrame(orf_results)
            
            # ORF Haritası
            st.pyplot(plot_orf_map(orf_results, len(dna_sequence)))
            
            # Sonuç Tablosu
            st.dataframe(df_orfs)
            
            # --- İNDİRME SEÇENEKLERİ ---
            st.subheader("📥 Analiz Sonuçlarını Dışa Aktar")
            col_dl1, col_dl2 = st.columns(2)
            
            # CSV İndirme
            csv = df_orfs.to_csv(index=False).encode('utf-8')
            col_dl1.download_button(label="CSV Olarak İndir", data=csv, file_name='orf_sonuclari.csv', mime='text/csv')
            
            
            # Excel İndirme
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df_orfs.to_excel(writer, index=False, sheet_name='ORF_Analizi')
            col_dl2.download_button(label="Excel (.xlsx) Olarak İndir", data=buffer.getvalue(), file_name='orf_sonuclari.xlsx', mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')

            st.divider()

            # 4. BİYOLOJİK ANOTASYON (BLAST)
            st.header("4. Protein Tanımlama ve Anotasyon (NCBI BLASTp)")
            orf_options = [f"ORF {i+1} ({row['AA Uzunluk']} aa, Çerçeve: {row['Çerçeve']})" for i, row in df_orfs.iterrows()]
            selected_orf_idx = st.selectbox("Anotasyon yapılacak ORF'yi seçin:", range(len(orf_options)), format_func=lambda x: orf_options[x])
            
            if st.button("Seçili ORF İçin BLAST Çalıştır (Zaman alabilir)"):
                selected_protein = df_orfs.iloc[selected_orf_idx]['Protein Dizisi']
                
                with st.spinner("NCBI veritabanı taranıyor..."):
                    try:
                        result_handle = NCBIWWW.qblast("blastp", "nr", selected_protein, hitlist_size=1)
                        blast_record = NCBIXML.read(result_handle)
                        
                        if blast_record.alignments:
                            alignment = blast_record.alignments[0]
                            hsp = alignment.hsps[0]
                            
                            st.success("Anotasyon Başarılı!")
                            st.markdown(f"**Protein/Gen Tanımı:** {alignment.title}")
                            st.markdown(f"**E-değeri (Güvenilirlik):** {hsp.expect}")
                            st.markdown(f"**Uzunluk:** {alignment.length} aa")
                            
                            with st.expander("Detaylı Hizalama (Alignment)"):
                                st.text(f"Query: {hsp.query}\nMatch: {hsp.match}\nSubjc: {hsp.sbjct}")
                        else:
                            st.warning("Anlamlı bir eşleşme bulunamadı.")
                    except Exception as e:
                        st.error(f"BLAST hatası: {e}")
                        
        else:
            st.warning("Belirtilen kriterlere uygun ORF bulunamadı.")
            
    else:
        st.error("Hata: Geçersiz DNA dizisi!")
else:
    st.info("Lütfen sol panelden bir veri girişi yapın.")