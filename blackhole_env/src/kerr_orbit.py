from __future__ import annotations
from dataclasses import dataclass
import math


'''
Aqui temos a representação da parte geométrica inicial do modelo. 
Para um buraco negro de Kerr com spin a, calculamos o horizonte externo e a ISCO
corrotante. Em seguida, verificamos se o raio orbital escolhido para o planeta está
fora do horizonte e fora da região onde órbitas circulares deixam de ser estáveis

G é a constante gravitacional, 
M é a massa do buraco negro,
c é a velocidade da luz

r+ é o horizonte externo,
r_isco é a ISCO corrotante,
z1 e z2 são variáveis auxiliares para o cálculo do risco de estabilidade
(contas intermediárias pro calculo da isco não ficar gigante)

a órbita é medida em tamanhos gravitacionais do buraco negro
(por isso adimensional em unidades geométricas G=c=M=1)

Como M=1 e a massa do buraconegro define a escala do sistema, r_orb é dado em função
de GM/c², que vai ser 1
'''


@dataclass(frozen=True)
class KerrOrbit:
    spin: float   # a, adimensional, em unidades geométricas G=c=M=1
    r_orb: float  # coordenada radial orbital, em unidades GM/c²

    # um Kerr com a>1 não é o caso físico de buraco negro de Kerr com horizonte
    def __post_init__(self):
        if not (0.0 <= self.spin <= 1.0):
            raise ValueError(
                f"Spin inválido: a={self.spin}. Para Kerr físico, use 0 <= a <= 1."
                )
        
        if self.r_orb <= 0:
            raise ValueError(
                f"Raio orbital inválido: r_orb={self.r_orb}. Deve ser positivo."
                )


    @property
    def r_plus(self) -> float:
        """Horizonte externo R+ = 1 + sqrt(1-a^2)."""
        a = self.spin
        return 1.0 + math.sqrt(max(0.0, 1.0 - a * a))

    @property
    def z1(self) -> float:
        a = self.spin
        return 1.0 + (1.0 - a * a) ** (1.0 / 3.0) * (
            (1.0 + a) ** (1.0 / 3.0) + (1.0 - a) ** (1.0 / 3.0)
        )

    @property
    def z2(self) -> float:
        return math.sqrt(3.0 * self.spin * self.spin + self.z1 * self.z1)

    def risco_corotating(self) -> float:
        """ISCO corrotante """
        z1 = self.z1
        z2 = self.z2
        a = self.spin
        return 3.0 + z2 - math.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))

# Verifica se a órbita é estável.
# Caso o raio da órbita não seja maior que o raio da ISCO corrotante (passou do limite de estabilidade)
#, lança um erro. (É necessário para que o modelo seja estável)
    def check_stable(self) -> None:
        risco = self.risco_corotating()
        if self.r_orb < risco:
            raise ValueError(
                f"Órbita dentro da ISCO corrotante: r_orb={self.r_orb:.8f}, "
                f"r_isco={risco:.8f}"
            )