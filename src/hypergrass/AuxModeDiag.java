package hypergrass;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

/** Reports the auxiliary coder's chosen mode on real channel dumps. Diagnostic only. */
public final class AuxModeDiag {
    public static void main(String[] args) throws Exception {
        Path dir = Path.of(args[0]);
        String[] names = args.length > 1
                ? Arrays.copyOfRange(args, 1, args.length)
                : new String[] {"gap", "diff", "mask", "geom", "litLen", "lits"};
        for (String name : names) {
            Path p = dir.resolve(name + ".bin");
            if (!Files.exists(p)) {
                System.out.println(name + ": missing");
                continue;
            }
            byte[] raw = Files.readAllBytes(p);
            VarIntCodec.ByteOut out = new VarIntCodec.ByteOut();
            long t0 = System.nanoTime();
            VarIntCodec.writeCompressedAux(out, raw, raw.length);
            long ms = (System.nanoTime() - t0) / 1_000_000;
            byte[] enc = Arrays.copyOf(out.internalBuffer(), out.size());
            byte[] back = VarIntCodec.readCompressedAux(new VarIntCodec.ByteIn(enc));
            System.out.println(name + ": raw=" + raw.length + " enc=" + enc.length
                    + " mode=" + (enc[0] & 0xFF) + " roundtrip=" + Arrays.equals(raw, back)
                    + " ms=" + ms);
        }
    }
}
