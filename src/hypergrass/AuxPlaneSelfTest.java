package hypergrass;

import java.util.Random;

/**
 * Round-trip test for the auxiliary-channel codecs, including the LEB128 digit-plane transform.
 * Verifies byte-exact recovery for varint streams, near-varint streams with a dirty tail, and
 * incompressible random data.
 */
public final class AuxPlaneSelfTest {

    public static void main(String[] args) {
        Random rnd = new Random(12345);
        check("tiny", new byte[] {1, 2, 3});
        check("empty", new byte[0]);

        for (int trial = 0; trial < 12; trial++) {
            int n = 2000 + rnd.nextInt(40000);
            byte[] raw = new byte[n];
            int pos = 0;
            int mode = trial % 4;
            while (pos < n) {
                int v;
                switch (mode) {
                    case 0 -> v = rnd.nextInt(4);
                    case 1 -> v = rnd.nextInt(64);
                    case 2 -> v = (int) Math.min(Integer.MAX_VALUE, Math.abs(rnd.nextGaussian()) * 4096);
                    default -> v = rnd.nextInt(1 << 20);
                }
                while ((v & ~0x7F) != 0) {
                    if (pos >= n) {
                        break;
                    }
                    raw[pos++] = (byte) ((v & 0x7F) | 0x80);
                    v >>>= 7;
                }
                if (pos < n) {
                    raw[pos++] = (byte) v;
                }
            }
            byte[] exact = java.util.Arrays.copyOf(raw, pos);
            check("varint-" + trial, exact);
        }

        byte[] noise = new byte[30000];
        rnd.nextBytes(noise);
        check("noise", noise);

        byte[] dirty = new byte[5000];
        rnd.nextBytes(dirty);
        check("dirty", dirty);

        // Long-range repetition: the pattern must reappear beyond Deflate's 32 KB window.
        byte[] block = new byte[900];
        for (int i = 0; i < block.length; i++) {
            block[i] = (byte) (i * 37 + 11);
        }
        byte[] far = new byte[200000];
        for (int i = 0; i < far.length; i++) {
            far[i] = block[i % block.length];
        }
        check("repeat-far", far);

        // Cross-block dictionary: a chunk that is entirely a copy of the prior 256 KB must
        // round-trip through the history-aware LZ coder and compress to a tiny token stream.
        byte[] dict = new byte[300000];
        for (int i = 0; i < dict.length; i++) {
            dict[i] = (byte) ((i * 31 + (i >>> 7)) & 0xFF);
        }
        byte[] hist = java.util.Arrays.copyOfRange(dict, 0, 256 * 1024);
        byte[] chunk = java.util.Arrays.copyOfRange(dict, 200000, 240000);
        byte[] dictEnc = LzRangeCoder.encode(chunk, chunk.length, hist, hist.length);
        byte[] dictBack = LzRangeCoder.decode(dictEnc, chunk.length, hist, hist.length);
        if (!java.util.Arrays.equals(chunk, dictBack)) {
            throw new IllegalStateException("history round-trip mismatch");
        }
        if (dictEnc.length * 10 >= chunk.length) {
            throw new IllegalStateException("history dictionary did not help: "
                    + chunk.length + " -> " + dictEnc.length);
        }
        System.out.println("history-dict: " + chunk.length + " -> " + dictEnc.length);

        byte[] mixed = new byte[120000];
        for (int i = 0; i < mixed.length; i++) {
            mixed[i] = (byte) (rnd.nextInt(100) < 70 ? block[i % block.length] : rnd.nextInt(256));
        }
        check("repeat-mixed", mixed);

        System.out.println("AuxPlaneSelfTest passed");
    }

    private static void check(String name, byte[] raw) {
        VarIntCodec.ByteOut out = new VarIntCodec.ByteOut();
        VarIntCodec.writeCompressedAux(out, raw, raw.length);
        byte[] enc = java.util.Arrays.copyOf(out.internalBuffer(), out.size());
        byte[] back = VarIntCodec.readCompressedAux(new VarIntCodec.ByteIn(enc));
        if (!java.util.Arrays.equals(raw, back)) {
            throw new IllegalStateException("round-trip mismatch: " + name
                    + " raw=" + raw.length + " enc=" + enc.length + " back=" + back.length);
        }
        System.out.println(name + ": " + raw.length + " -> " + enc.length + " mode=" + (enc[0] & 0xFF));
    }
}
