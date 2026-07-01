# Resumo — Alterações e Conclusões

**Romeu Fernandes (PG45861) · Universidade do Minho · 2026-06-16**

---

## O que alterei

- **Conservação → PPI** (Percentage of Pairwise Identity, do ViruScope / Ana Lima):
  identidade par-a-par IUPAC-aware, em vez da simples "base mais comum".
- **No-fold score** (ViruScope): converte o ΔG num score contínuo 0–100, em vez de só pass/fail.
- **Seleção do comprimento automática por gene** (cluster dominante pós-fetch): resolveu o
  **algD (0 → 403 probes)** e dispensa afinar `min_len`/`max_len` à mão em genes novos.
- **Limiares recalibrados com base na literatura e nos dados:** `gc_max 0.65 → 0.60`;
  **`seqfold −6,0 → −2,6 kcal/mol` (= percentil 5 da distribuição observada)**.
- **Janela de Tm automática** a partir da temperatura do ensaio (`--assay-temp`).
- **Nº de sequências (N=100) justificado por rarefação** (saturação da diversidade).
- **Seleção das probes para o Boltz por qualidade** (PPI + No-fold), já não pelo Tm.
- **Análise exploratória** completa: distribuição do seqfold, tamanhos + normalidade,
  diversidade sem alinhamento, rarefação e descritores de sequência.
- **Estrutura 3D (Boltz-2)** com **matriz PAE** (estilo AlphaFold3) + shortlist final ranqueada.
- **Probes externas (IPLEX)** integradas pelo fluxo normal, pontuadas com as mesmas métricas.

## Conclusões

- **Funil de triagem:** 8 455 candidatas → **3 078 passam a triagem básica** →
  **2 943 passam o seqfold** → **30 selecionadas** para previsão 3D.
- **Validação 3D (Boltz):** **4 probes GOOD**, 11 MODERATE, 15 LOW. As melhores são sólidas
  **em sequência E em 3D**:

  | probe | gene | PPI | No-fold | confidence | pLDDT |
  |---|---|---|---|---|---|
  | p5972 | algD | 0,93 | 98 | **0,731** | 0,901 |
  | p7716 | frdB | 0,99 | 100 | 0,723 | 0,822 |
  | p5938 | algD | 0,90 | 96 | 0,721 | 0,887 |
  | p8448 | frdB | 1,00 | 98 | 0,716 | 0,845 |

- O **limiar de seqfold −2,6 kcal/mol é o percentil 5** da distribuição real das probes →
  remove os ~5% com estrutura secundária mais estável (decisão fundamentada, não arbitrária).
- **N ≈ 100 sequências/gene é suficiente** (a diversidade satura por volta de 50–75); o fator
  limitante real é o número de sequências no NCBI (ex.: algD, oprL), não o limite escolhido.
- O **algD**, que antes não dava probes, passou a dar 403 e produziu **as duas melhores
  estruturas 3D** — valida a seleção automática de comprimento.
- **pTM baixo é esperado** para ssDNA curto (a métrica é calibrada para proteínas); os
  indicadores úteis aqui são o *confidence* e o *pLDDT*.

*Detalhe e referências em [`docs/Justificacao_Pipeline.md`](Justificacao_Pipeline.md); dados e figuras em `docs/analysis/`.*
