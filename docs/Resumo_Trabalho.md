# Resumo — Alterações e conclusões

**Romeu Fernandes (PG45861) · Universidade do Minho**

---

## 1. Parâmetros por tipo de organismo e por espécie

**Porquê (a lógica).** O painel de alvos não é só bactérias — inclui **vírus, fungos,
protozoários e um controlo humano**. A composição do genoma (**GC**) e a **variabilidade**
diferem enormemente entre eles, por isso os mesmos limiares não servem para todos:
- ***Plasmodium falciparum*** (protozoário, malária): genoma **~19% GC** (AT-rico) → as suas
  probes são AT-ricas; os critérios bacterianos (GC 40–60%) **rejeitá-las-iam injustamente**.
- **Herpes (HSV-1/2):** ~68–70% GC (GC-rico); **influenza (H1N1)** ~43% GC → os vírus precisam
  de uma **janela GC larga**.
- **Fungo** (*Aspergillus*) ~50% GC; **hospedeiro humano** ~41% GC.

**O que fiz.** Criei **perfis de parâmetros por tipo** (bactéria / vírus / fungo / protozoário /
hospedeiro) e **por espécie** (afinados pelo GC do genoma de cada uma). A resolução é em camadas:
**gene → espécie → tipo → global** (o mais específico ganha). Os vírus, por serem mais variáveis,
têm ainda uma **conservação mais permissiva**.

**Espécies novas — definidas pelo utilizador em tempo real.** Quando o pipeline encontra uma
espécie sem perfil, **pergunta** o tipo e cada limiar, **sugerindo** um valor (que se aceita com
Enter ou se altera), e guarda em `data/species_params.yaml` para as próximas vezes.

**Resultado.** Ao pontuar as 74 probes IPLEX pelos critérios da **sua** espécie/tipo,
**29 → 39** passam a triagem básica — ex.: a probe de *Plasmodium* (GC 0,30, AT-rica) **passa
agora**, o que os critérios bacterianos globais rejeitavam. As referências de cada perfil estão em
`docs/parametros_referencias.csv` (para a secção de Métodos da tese).

## 2. As probes IPLEX no pipeline e a comparação 3D

- As 74 probes de `data/sequencias_iplex.xlsx` entram pelo **fluxo normal** (`--with-reference`),
  pontuadas com as **mesmas métricas** das minhas (Tm, GC, hairpin, homodímero, seqfold, No-fold)
  e pelos critérios da **sua espécie** (§1). Ficam no `FINAL_PROBES_ALL.csv` com `fonte = referencia`,
  **lado a lado com as minhas** e com as mesmas colunas → comparação direta na mesma tabela.
- A **conservação (PPI) fica não aplicável** — só é definível para probes desenhadas a partir de um
  alinhamento; probes já feitas não têm região de origem.
- Em média, as IPLEX que passam ficam **ligeiramente abaixo** das minhas (Tm ~56 vs 59; No-fold
  ~80 vs 83) — coerente com serem um painel mais antigo e heterogéneo.
- **Boltz 3D (um comando):** `python pipeline.py --max-seqs 100 --with-reference --colab 5` gera o
  input com **as minhas top-5/gene + as IPLEX que passam** (68 probes = 30 + 38) → previsão 3D no
  Colab → `--merge-boltz` junta tudo numa shortlist para comparar **minha vs IPLEX em 3D**.

## 3. Novas métricas de diversidade (sem alinhamento)

Todas são **baseadas em k-mers, sem alinhamento** (o objetivo era avaliar a diversidade *sem
alinhar* as sequências), pelo que funcionam mesmo com comprimentos variáveis e sequências divergentes.

### 3.1 Porquê estas métricas (e não outras)
São **complementares** — cada uma mede um aspeto diferente e, em conjunto, dão uma leitura robusta
que uma média sozinha não daria:
- **Cosseno** pesa as *frequências* dos k-mers; **Jaccard** usa só *presença/ausência* (o
  "vocabulário" de k-mers). Se as duas concordam, a diversidade é real, não um artefacto de uma métrica.
- **Desvio-padrão** e **máximo** do cosseno mostram a *forma* da distribuição (subgrupos/outliers
  que a média esconde).
- **% únicas** deteta o artefacto dos *duplicados* (o NCBI devolve muitas cópias).

Preteri a distância **Euclidiana** (sensível ao comprimento das sequências) e as medidas
**baseadas em alinhamento** (π, identidade média — violam o requisito "sem alinhamento", são lentas
e falham em conjuntos divergentes).

### 3.2 Como ler cada métrica

| Métrica | O que mede | Valor **alto** → | Valor **baixo** → |
|---|---|---|---|
| Distância de cosseno média | distância composicional média (frequências de k-mers) | conjunto **diverso** | conjunto **conservado/semelhante** |
| Desvio-padrão da distância de cosseno | dispersão da diversidade entre pares | **heterogéneo** (há subgrupos/outliers) | diversidade **uniforme** |
| Distância de cosseno máxima | o par de sequências mais divergente | existe pelo menos um **outlier**/subgrupo distinto | até o par mais diferente é **próximo** |
| Distância de Jaccard média | distância por presença/ausência de k-mers | partilham **poucos** k-mers (vocabulário diferente) | partilham **quase todos** os k-mers |
| Sequências únicas (%) | redundância (duplicados exatos) | sequências **todas distintas** | **muitos duplicados** no conjunto |

### 3.3 Resultados (por gene)

| gene | nº seqs | cosseno médio | desvio-padrão | cosseno máximo | Jaccard médio | sequências únicas (%) |
|---|---|---|---|---|---|---|
| nuc  | 77 | 0,059 | 0,074 | 0,425 | 0,163 | 61,0 |
| rmpM | 76 | 0,079 | 0,076 | 0,247 | 0,069 | 50,0 |
| lytA | 87 | **0,187** | 0,144 | 0,440 | 0,125 | 93,1 |
| oprL | 24 | 0,050 | 0,084 | 0,358 | 0,096 | 91,7 |
| algD | 12 | 0,097 | 0,140 | 0,332 | 0,117 | **100,0** |
| frdB | 82 | 0,025 | 0,013 | 0,052 | 0,061 | **35,4** |

**Conclusões:**
- **lytA** é o mais diverso (cosseno médio 0,187) — coerente com a sua diversidade alélica.
- **frdB** é o mais conservado e redundante (cosseno 0,025 e só **35% de sequências únicas** →
  muitos duplicados no NCBI).
- **algD** tem 100% de sequências únicas (poucas, mas todas distintas) → explica a sua instabilidade
  e reforça a necessidade da seleção automática de comprimento.
- **Cosseno e Jaccard concordam** entre si → leitura robusta da diversidade, sem alinhamento.
