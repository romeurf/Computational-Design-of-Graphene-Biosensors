# O que foi implementado desde a última reunião

**Romeu Fernandes (PG45861) · Universidade do Minho** · 2026-07-13

Na reunião ficaram quatro pedidos: (1) parâmetros diferentes para vírus e bactérias;
(2) valores otimizados por espécie; (3) um perfil geral, mas deixando o **utilizador definir
os limites** para espécies novas; (4) correr as probes do **IPLEX** juntamente com as minhas no
**Boltz** (pipeline completa). Está tudo implementado no `pipeline.py` e validado numa corrida real.
Segue o que cada ponto significa e como ficou.

---

## 1. Parâmetros por tipo de organismo e por espécie

**Porquê.** O painel de alvos não é só bactérias — inclui **vírus, fungos, protozoários e um
controlo humano**. A composição do genoma (**%GC**) e a **variabilidade** diferem enormemente,
por isso os mesmos limiares não servem para todos. Exemplos concretos:
- ***Plasmodium falciparum*** (protozoário, malária): genoma **~19% GC** (AT-rico) → as suas probes
  são AT-ricas; os critérios bacterianos (GC 40–60%) **rejeitá-las-iam injustamente**.
- **Herpes (HSV-1/2):** ~68–70% GC (GC-rico); **influenza (H1N1)** ~43% GC → os vírus precisam de
  uma **janela GC larga**.
- **Fungo** (*Aspergillus*) ~50% GC; **hospedeiro humano** ~41% GC.

**O que fiz.** Criei **perfis de parâmetros em três camadas**, do mais específico ao mais geral:

> **gene → espécie → tipo → global** *(o mais específico ganha)*

- **Por tipo** (`TYPE_DEFAULTS`): 5 perfis — **bactéria / vírus / fungo / protozoário / hospedeiro**.
  As diferenças vêm sobretudo do GC do genoma e da variabilidade; os vírus, por serem mais variáveis,
  têm ainda uma **conservação mais permissiva**.
- **Por espécie** (`SPECIES_PARAMS`): 17 espécies afinadas pelo GC do seu genoma (HSV-1 ~68%,
  *P. aeruginosa* ~67%, *Plasmodium* ~19%, etc.).

Cada perfil guarda as **referências bibliográficas** que justificam os valores — reunidas num ficheiro
consultável, `docs/parametros_referencias.csv`, para a secção de Métodos da tese.

---

## 2. Espécies novas: definidas pelo utilizador em tempo real

Este é o ponto dos **inputs novos**. Em vez de uma espécie desconhecida cair em silêncio num default
qualquer, o pipeline **pára e pergunta no momento**, sugerindo sempre um valor de partida:

```
  → Organismo sem perfil: 'Influenza A H3N2'
    Tipo (bacteria/virus/fungus/protozoa/host)? [sugestão: virus]  ⏎
    GC mínimo [sugestão 0.35]: 0.30
    GC máximo [sugestão 0.68]: ⏎
    ...
    ✔ perfil de 'Influenza A H3N2' guardado em species_params.yaml
```

- **Primeiro define-se o tipo**, depois cada limiar. O programa **sugere** um valor (a partir do tipo e
  do GC do genoma); carregar **Enter aceita**, ou escreve-se outro valor.
- O tipo sugerido vem de um heurístico (`_infer_type`): procura a espécie no catálogo e, se não a
  encontrar, deduz pelo nome (p. ex. "*Human herpesvirus*" → vírus, "*Aspergillus*" → fungo).
- O perfil resultante é **guardado em `data/species_params.yaml`** → nas próximas corridas já não volta
  a perguntar por essa espécie.
- **Não bloqueia em modo automático:** só pergunta se houver terminal interativo; numa corrida em lote
  usa as sugestões silenciosamente.

---

## 3. Juntar as probes do IPLEX com as minhas

As **74 probes do IPLEX** (`data/sequencias_iplex.xlsx`, 17 espécies) entram agora pelo **fluxo normal**
do pipeline (opção `--with-reference`):
- São pontuadas com as **mesmas métricas** das minhas (Tm, GC, hairpin, homodímero, seqfold, No-fold) e
  **pelos critérios da sua própria espécie/tipo** (§1) — não pelos critérios bacterianos globais.
- Ficam no `FINAL_PROBES_ALL.csv` marcadas com `fonte = referencia`, **lado a lado com as minhas** e com
  as mesmas colunas → comparação direta na mesma tabela.
- A **conservação fica "não aplicável"** para elas — só é definível para probes desenhadas a partir de um
  alinhamento; probes já feitas não têm região de origem.

**Efeito dos critérios por espécie:** ao pontuar as IPLEX pelos limiares da *sua* espécie, mais probes
passam a triagem básica do que passariam com os critérios bacterianos globais — por exemplo, a probe de
*Plasmodium* (AT-rica) **passa agora**, quando antes seria rejeitada.

---

## 4. Correr IPLEX + as minhas no Boltz (pipeline completa)

Um único comando prepara a comparação 3D com **as duas fontes juntas**:

```bash
python pipeline.py --max-seqs 100 --with-reference --colab 5
```

Gera um `boltz2_inputs.zip` com **as minhas top-5 por gene + as IPLEX que passam o crivo**, cada uma
identificada pela `fonte`. Sobe-se ao Colab (GPU) para a previsão 3D e depois:

```bash
python pipeline.py --merge-boltz output/colab_boltz2/boltz2_results_summary.csv
```

junta os resultados 3D numa shortlist (minha vs IPLEX no mesmo ranking) **e** regenera automaticamente o
ficheiro de entrega para a pipeline de docking externa (`docs/melhores_probes_docking.csv`, no formato
`job_name, sequence, species`), contendo **só as minhas melhores probes**, ordenadas por confiança 3D.

---

## Validação — corrida real (2026-07-13)

Corrida completa de ponta a ponta, sem erros:

| Etapa | Resultado |
|---|---|
| Genes recuperados do NCBI | 6/6 (nuc, rmpM, lytA, oprL, algD, frdB) |
| Funil global | 8455 janelas → **3078** passam triagem básica → **2943** passam seqfold |
| Probes IPLEX pontuadas | **74** (pelos critérios da sua espécie) → **38 passam** |
| Tabela consolidada | `FINAL_PROBES_ALL.csv` — **8529 probes** (8455 minhas + 74 IPLEX) |
| Input Boltz | `boltz2_inputs.zip` — **68 probes** (30 minhas top-5/gene + 38 IPLEX) |

As 17 espécies do painel IPLEX estão todas pré-definidas, por isso esta corrida não precisou de perguntar
nada; o mecanismo de definição interativa (§2) só é acionado por espécies **novas** fora do catálogo.

---

## Ficheiros relevantes

- `pipeline.py` — perfis por tipo/espécie, resolução em camadas, definição interativa, IPLEX no fluxo,
  export para Boltz, entrega para docking.
- `data/species_params.yaml` — perfis de espécies definidos pelo utilizador (criado/editável).
- `docs/parametros_referencias.csv` — referências bibliográficas de cada perfil (para os Métodos).
- `docs/melhores_probes_docking.csv` — entrega para a pipeline de docking (auto-atualizada).
